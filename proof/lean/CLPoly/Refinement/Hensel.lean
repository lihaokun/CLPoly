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
import CLPoly.Refinement.DDF
import CLPoly.Refinement.PolynomialMod
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
  toPolyMod m output.t = toPolyMod m input.t ∧
  StrictPolynomialMod.SparsePolyZZCanonical output.g ∧
  StrictPolynomialMod.SparsePolyZZCanonical output.h ∧
  ∀ head tail, input.h.toList = (head, 1) :: tail → 1 < m →
    ∃ suffix, output.h.toList = (head, 1) :: suffix

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

/-- Signed-integer view of the same exact source degree read. -/
theorem get_deg_toInt_eq_head (f : SparsePolyZZ)
    (hnonempty : 0 < f.size) (hdegree : f[0]!.1.deg < 2 ^ 63) :
    (get_deg f).toInt = (f[0]!.1.deg : Int) := by
  have hne : f ≠ #[] := by
    intro hempty
    subst f
    simp at hnonempty
  have hget : get_deg f = f[0]!.1.deg.toUInt64.toInt64 := by
    simp [get_deg, Array.isEmpty_iff, hne, getElem!_pos f 0 hnonempty]
  rw [hget]
  change f[0]!.1.deg.toInt64.toInt = (f[0]!.1.deg : Int)
  exact Int64.toInt_ofNat_of_lt hdegree

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

/-- Coefficient reduction only filters terms and leaves every retained
degree untouched, so it preserves the source `upolynomial_` representation
invariant needed by modular long division. -/
theorem modCoeffOutput_canonical (f : SparsePolyZZ) (m : Int)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.modCoeffOutput f m) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  constructor
  · let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
      fun a b => a.1.deg > b.1.deg
    letI : Trans degreeGreater degreeGreater degreeGreater :=
      ⟨by
        intro a b c hab hbc
        dsimp [degreeGreater] at hab hbc ⊢
        omega⟩
    apply List.isChain_iff_pairwise.mpr
    rw [Generated.StrictHensel.modCoeffOutput, Array.toList_filterMap]
    apply List.Pairwise.filterMap (R := degreeGreater) (S := degreeGreater)
      (fun term =>
        let coefficient := ZZ.fdiv_r term.2 term.2 m
        if coefficient != 0 then some (term.1, coefficient) else none)
      (fun a b hab outputA ha outputB hb => by
        dsimp only at ha hb
        split at ha <;> try contradiction
        split at hb <;> try contradiction
        simp only [Option.some.injEq] at ha hb
        subst outputA
        subst outputB
        exact hab)
    exact List.isChain_iff_pairwise.mp hf.1
  · intro term hterm
    rw [Generated.StrictHensel.modCoeffOutput,
      Array.toList_filterMap, List.mem_filterMap] at hterm
    rcases hterm with ⟨source, _, hsource⟩
    dsimp only at hsource
    split at hsource
    next hnonzero =>
      simp only [Option.some.injEq] at hsource
      subst term
      simpa using hnonzero
    next hzero => contradiction

/-- The generated coefficient-scaling range-for preserves the canonical
sparse representation when its integer scalar is nonzero. -/
theorem scaleCoeffs_canonical (f : SparsePolyZZ) (scale : Int)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hscale : scale ≠ 0) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.scaleCoeffs f scale) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  constructor
  · rw [Generated.StrictHensel.scaleCoeffs, Array.toList_map,
      List.isChain_map]
    simpa using hf.1
  · intro term hterm
    rw [Generated.StrictHensel.scaleCoeffs, Array.toList_map,
      List.mem_map] at hterm
    rcases hterm with ⟨source, hsource, rfl⟩
    exact mul_ne_zero (hf.2 source hsource) hscale

/-- The generated divide/reduce coefficient loop preserves strict degree
order and removes every zero coefficient it creates. -/
theorem divideThenReduceCoeffs_canonical (f : SparsePolyZZ) (modulus : Int)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.divideThenReduceCoeffs f modulus) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  constructor
  · let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
      fun a b => a.1.deg > b.1.deg
    letI : Trans degreeGreater degreeGreater degreeGreater :=
      ⟨by
        intro a b c hab hbc
        dsimp [degreeGreater] at hab hbc ⊢
        omega⟩
    apply List.isChain_iff_pairwise.mpr
    rw [Generated.StrictHensel.divideThenReduceCoeffs,
      Array.toList_filterMap]
    apply List.Pairwise.filterMap (R := degreeGreater) (S := degreeGreater)
      (fun term =>
        let quotient := ZZ.fdiv_q term.2 term.2 modulus
        let coefficient := ZZ.fdiv_r quotient quotient modulus
        if coefficient != 0 then some (term.1, coefficient) else none)
      (fun a b hab outputA ha outputB hb => by
        dsimp only at ha hb
        split at ha <;> try contradiction
        split at hb <;> try contradiction
        simp only [Option.some.injEq] at ha hb
        subst outputA
        subst outputB
        exact hab)
    exact List.isChain_iff_pairwise.mp hf.1
  · intro term hterm
    rw [Generated.StrictHensel.divideThenReduceCoeffs,
      Array.toList_filterMap, List.mem_filterMap] at hterm
    rcases hterm with ⟨source, _, hsource⟩
    dsimp only at hsource
    split at hsource
    next hnonzero =>
      simp only [Option.some.injEq] at hsource
      subst term
      simpa using hnonzero
    next hzero => contradiction

/-- The exact composition used by the final Hensel normalization is canonical
whenever its source factor is canonical.  This is proved directly for the
generated map/filter traversal: the intermediate scaled array may contain
zero coefficients and therefore need not itself satisfy the invariant. -/
theorem modCoeffOutput_scaleCoeffs_canonical (f : SparsePolyZZ)
    (scale modulus : Int)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.modCoeffOutput
        (Generated.StrictHensel.scaleCoeffs f scale) modulus) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  constructor
  · let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
      fun a b => a.1.deg > b.1.deg
    letI : Trans degreeGreater degreeGreater degreeGreater :=
      ⟨by
        intro a b c hab hbc
        dsimp [degreeGreater] at hab hbc ⊢
        omega⟩
    apply List.isChain_iff_pairwise.mpr
    rw [Generated.StrictHensel.modCoeffOutput, Array.toList_filterMap,
      Generated.StrictHensel.scaleCoeffs, Array.toList_map]
    apply List.Pairwise.filterMap (R := degreeGreater) (S := degreeGreater)
      (fun term =>
        let coefficient := ZZ.fdiv_r term.2 term.2 modulus
        if coefficient != 0 then some (term.1, coefficient) else none)
      (fun a b hab outputA ha outputB hb => by
        dsimp only at ha hb
        split at ha <;> try contradiction
        split at hb <;> try contradiction
        simp only [Option.some.injEq] at ha hb
        subst outputA
        subst outputB
        exact hab)
    rw [List.pairwise_map]
    simpa [degreeGreater] using List.isChain_iff_pairwise.mp hf.1
  · intro term hterm
    rw [Generated.StrictHensel.modCoeffOutput, Array.toList_filterMap,
      Generated.StrictHensel.scaleCoeffs, Array.toList_map,
      List.mem_filterMap] at hterm
    rcases hterm with ⟨source, _, hsource⟩
    dsimp only at hsource
    split at hsource
    next hnonzero =>
      simp only [Option.some.injEq] at hsource
      subst term
      simpa using hnonzero
    next hzero => contradiction

/-- Removing the leading array cell preserves strict degree order and all
remaining nonzero coefficients. -/
theorem eraseLeading_canonical (f : SparsePolyZZ) (hfsize : 0 < f.size)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    StrictPolynomialMod.SparsePolyZZCanonical (f.eraseIdxIfInBounds 0) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical at hf ⊢
  rw [Array.toList_eraseIdxIfInBounds, List.eraseIdx_zero]
  exact ⟨hf.1.tail, fun term hterm => hf.2 term (List.mem_of_mem_tail hterm)⟩

/-- In a nonempty canonical sparse polynomial, every stored term after the
front cell has strictly smaller degree than the front cell. -/
theorem canonical_tail_degree_lt_head (f : SparsePolyZZ) (hfsize : 0 < f.size)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    ∀ term ∈ f.toList.drop 1, term.1.deg < f[0]!.1.deg := by
  have hfull : f.toList = f[0]! :: f.toList.drop 1 := by
    rw [getElem!_pos f 0 hfsize]
    simpa using List.drop_eq_getElem_cons
      (l := f.toList) (i := 0) (by simpa using hfsize)
  have hpairwise : List.Pairwise
      (fun a b : UMonomial × Int => a.1.deg > b.1.deg) f.toList := by
    let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
      fun a b => a.1.deg > b.1.deg
    letI : Trans degreeGreater degreeGreater degreeGreater :=
      ⟨by
        intro a b c hab hbc
        dsimp [degreeGreater] at hab hbc ⊢
        omega⟩
    exact List.isChain_iff_pairwise.mp hf.1
  rw [hfull, List.pairwise_cons] at hpairwise
  exact hpairwise.1

theorem eraseLeading_size_lt (f : SparsePolyZZ) (hfsize : 0 < f.size) :
    (f.eraseIdxIfInBounds 0).size < f.size := by
  unfold Array.eraseIdxIfInBounds
  rw [dif_pos hfsize, Array.size_eraseIdx]
  omega

private theorem pushNonzero_canonical (result : SparsePolyZZ) (degree : Nat)
    (coefficient : Int)
    (hresult : StrictPolynomialMod.SparsePolyZZCanonical result)
    (habove : ∀ term ∈ result.toList, degree < term.1.deg) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.pushNonzero result degree coefficient) := by
  by_cases hzero : coefficient = 0
  · subst coefficient
    simpa [Generated.StrictHensel.pushNonzero] using hresult
  · unfold Generated.StrictHensel.pushNonzero
    simp only [bne_iff_ne]
    rw [if_pos hzero]
    unfold StrictPolynomialMod.SparsePolyZZCanonical at hresult ⊢
    rw [Array.toList_push]
    constructor
    · let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
        fun a b => a.1.deg > b.1.deg
      letI : Trans degreeGreater degreeGreater degreeGreater :=
        ⟨by
          intro a b c hab hbc
          dsimp [degreeGreater] at hab hbc ⊢
          omega⟩
      apply List.isChain_iff_pairwise.mpr
      rw [List.pairwise_append]
      refine ⟨List.isChain_iff_pairwise.mp hresult.1, by simp, ?_⟩
      intro term hterm added hadd
      simp only [List.mem_singleton] at hadd
      subst added
      exact habove term hterm
    · intro term hterm
      simp only [List.mem_append, List.mem_singleton] at hterm
      rcases hterm with hold | rfl
      · exact hresult.2 term hold
      · exact hzero

private theorem canonical_pairwise_drop (f : SparsePolyZZ) (index : Nat)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    (f.toList.drop index).Pairwise
      (fun a b : UMonomial × Int => b.1.deg < a.1.deg) := by
  let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
    fun a b => a.1.deg > b.1.deg
  letI : Trans degreeGreater degreeGreater degreeGreater :=
    ⟨by
      intro a b c hab hbc
      dsimp [degreeGreater] at hab hbc ⊢
      omega⟩
  exact (List.isChain_iff_pairwise.mp hf.1).drop

private theorem sparseDrop_eq_getElem_cons (f : SparsePolyZZ) (index : Nat)
    (hindex : index < f.size) :
    f.toList.drop index = f[index]! :: f.toList.drop (index + 1) := by
  have hlength : index < f.toList.length := by simpa using hindex
  have hdrop := List.drop_eq_getElem_cons
    (l := f.toList) (i := index) hlength
  rw [getElem!_pos f index hindex]
  simpa using hdrop

private theorem canonical_cursor_tail_lt (f : SparsePolyZZ) (index : Nat)
    (hindex : index < f.size)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    ∀ term ∈ f.toList.drop (index + 1),
      term.1.deg < f[index]!.1.deg := by
  have hordered := canonical_pairwise_drop f index hf
  rw [sparseDrop_eq_getElem_cons f index hindex,
    List.pairwise_cons] at hordered
  exact hordered.1

private theorem canonical_cursor_degree_le (f : SparsePolyZZ) (index : Nat)
    (hindex : index < f.size)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) :
    ∀ term ∈ f.toList.drop index, term.1.deg ≤ f[index]!.1.deg := by
  rw [sparseDrop_eq_getElem_cons f index hindex]
  intro term hterm
  simp only [List.mem_cons] at hterm
  rcases hterm with rfl | htail
  · exact le_rfl
  · exact (canonical_cursor_tail_lt f index hindex hf term htail).le

private structure PairVecAddCursorValid (a b : SparsePolyZZ)
    (aIndex bIndex : Nat) (result : SparsePolyZZ) : Prop where
  result_canonical : StrictPolynomialMod.SparsePolyZZCanonical result
  result_above_a : ∀ output ∈ result.toList,
    ∀ term ∈ a.toList.drop aIndex, term.1.deg < output.1.deg
  result_above_b : ∀ output ∈ result.toList,
    ∀ term ∈ b.toList.drop bIndex, term.1.deg < output.1.deg

private theorem PairVecAddCursorValid.pushNonzero
    (a b : SparsePolyZZ) (aIndex bIndex nextA nextB : Nat)
    (result : SparsePolyZZ) (degree : Nat) (coefficient : Int)
    (hvalid : PairVecAddCursorValid a b aIndex bIndex result)
    (hdegree : ∀ output ∈ result.toList, degree < output.1.deg)
    (haBelow : ∀ term ∈ a.toList.drop nextA, term.1.deg < degree)
    (hbBelow : ∀ term ∈ b.toList.drop nextB, term.1.deg < degree) :
    PairVecAddCursorValid a b nextA nextB
      (Generated.StrictHensel.pushNonzero result degree coefficient) := by
  refine ⟨pushNonzero_canonical result degree coefficient
      hvalid.result_canonical hdegree, ?_, ?_⟩
  · intro output houtput term hterm
    by_cases hzero : coefficient = 0
    · subst coefficient
      have hold : output ∈ result.toList := by
        simpa [Generated.StrictHensel.pushNonzero] using houtput
      exact lt_trans (haBelow term hterm) (hdegree output hold)
    · rw [Generated.StrictHensel.pushNonzero, if_pos (by simpa using hzero),
        Array.toList_push] at houtput
      simp only [List.mem_append, List.mem_singleton] at houtput
      rcases houtput with hold | rfl
      · exact lt_trans (haBelow term hterm) (hdegree output hold)
      · exact haBelow term hterm
  · intro output houtput term hterm
    by_cases hzero : coefficient = 0
    · subst coefficient
      have hold : output ∈ result.toList := by
        simpa [Generated.StrictHensel.pushNonzero] using houtput
      exact lt_trans (hbBelow term hterm) (hdegree output hold)
    · rw [Generated.StrictHensel.pushNonzero, if_pos (by simpa using hzero),
        Array.toList_push] at houtput
      simp only [List.mem_append, List.mem_singleton] at houtput
      rcases houtput with hold | rfl
      · exact lt_trans (hbBelow term hterm) (hdegree output hold)
      · exact hbBelow term hterm

private theorem PairVecAddCursorValid.push
    (a b : SparsePolyZZ) (aIndex bIndex nextA nextB : Nat)
    (result : SparsePolyZZ) (degree : Nat) (coefficient : Int)
    (hcoefficient : coefficient ≠ 0)
    (hvalid : PairVecAddCursorValid a b aIndex bIndex result)
    (hdegree : ∀ output ∈ result.toList, degree < output.1.deg)
    (haBelow : ∀ term ∈ a.toList.drop nextA, term.1.deg < degree)
    (hbBelow : ∀ term ∈ b.toList.drop nextB, term.1.deg < degree) :
    PairVecAddCursorValid a b nextA nextB
      (result.push (UMonomial.mk degree, coefficient)) := by
  simpa [Generated.StrictHensel.pushNonzero, hcoefficient] using
    hvalid.pushNonzero a b aIndex bIndex nextA nextB result degree coefficient
      hdegree haBelow hbBelow

set_option maxHeartbeats 0 in
/-- The exact generated `pair_vec_add` merge preserves the canonical sparse
integer representation at every cursor state satisfying the emitted-prefix
frontier invariant. -/
theorem pairVecAddLoop_canonical (a b : SparsePolyZZ)
    (ha : StrictPolynomialMod.SparsePolyZZCanonical a)
    (hb : StrictPolynomialMod.SparsePolyZZCanonical b) :
    ∀ aIndex bIndex result,
      PairVecAddCursorValid a b aIndex bIndex result →
      StrictPolynomialMod.SparsePolyZZCanonical
        (Generated.StrictHensel.pairVecAddLoop a b aIndex bIndex result) := by
  intro aIndex bIndex result hvalid
  refine Generated.StrictHensel.pairVecAddLoop.induct a b
    (motive := fun aIndex bIndex result =>
      PairVecAddCursorValid a b aIndex bIndex result →
      StrictPolynomialMod.SparsePolyZZCanonical
        (Generated.StrictHensel.pairVecAddLoop a b aIndex bIndex result))
    ?_ ?_ ?_ ?_ ?_ ?_ aIndex bIndex result hvalid
  · intro ai bi acc hmore haDone ih hcursor
    have hbi : bi < b.size := by omega
    have hbmem : b[bi]! ∈ b.toList := by
      rw [getElem!_pos b bi hbi]
      exact Array.getElem_mem_toList hbi
    have hcoeff : b[bi]!.2 ≠ 0 := hb.2 b[bi]! hbmem
    have haDrop : a.toList.drop ai = [] := List.drop_eq_nil_of_le (by
      simpa using haDone)
    have hnext := hcursor.push a b ai bi ai (bi + 1) acc b[bi]!.1.deg
      b[bi]!.2 hcoeff
      (fun output houtput => hcursor.result_above_b output houtput b[bi]!
        (by rw [drop_eq_getElem_cons b bi hbi]; simp))
      (by simp [haDrop])
      (canonical_cursor_tail_lt b bi hbi hb)
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp only [hmore, haDone, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hmore haMore hbDone ih hcursor
    have hai : ai < a.size := by omega
    have hamem : a[ai]! ∈ a.toList := by
      rw [getElem!_pos a ai hai]
      exact Array.getElem_mem_toList hai
    have hcoeff : a[ai]!.2 ≠ 0 := ha.2 a[ai]! hamem
    have hbDrop : b.toList.drop bi = [] := List.drop_eq_nil_of_le (by
      simpa using hbDone)
    have hnext := hcursor.push a b ai bi (ai + 1) bi acc a[ai]!.1.deg
      a[ai]!.2 hcoeff
      (fun output houtput => hcursor.result_above_a output houtput a[ai]!
        (by rw [drop_eq_getElem_cons a ai hai]; simp))
      (canonical_cursor_tail_lt a ai hai ha)
      (by simp [hbDrop])
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp only [hmore, haMore, hbDone, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hmore haMore hbMore hdegree ih hcursor
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have hbmem : b[bi]! ∈ b.toList := by
      rw [getElem!_pos b bi hbi]
      exact Array.getElem_mem_toList hbi
    have hcoeff : b[bi]!.2 ≠ 0 := hb.2 b[bi]! hbmem
    have hnext := hcursor.push a b ai bi ai (bi + 1) acc b[bi]!.1.deg
      b[bi]!.2 hcoeff
      (fun output houtput => hcursor.result_above_b output houtput b[bi]!
        (by rw [drop_eq_getElem_cons b bi hbi]; simp))
      (fun term hterm => lt_of_le_of_lt
        (canonical_cursor_degree_le a ai hai ha term hterm) hdegree)
      (canonical_cursor_tail_lt b bi hbi hb)
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp only [hmore, haMore, hbMore, hdegree, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hmore haMore hbMore hnotGreater hequal ih hcursor
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have habove : ∀ output ∈ acc.toList,
        a[ai]!.1.deg < output.1.deg := by
      intro output houtput
      exact hcursor.result_above_a output houtput a[ai]!
        (by rw [drop_eq_getElem_cons a ai hai]; simp)
    have hnext := hcursor.pushNonzero a b ai bi (ai + 1) (bi + 1) acc
      a[ai]!.1.deg (b[bi]!.2 + a[ai]!.2) habove
      (canonical_cursor_tail_lt a ai hai ha)
      (fun term hterm => by
        have hlt := canonical_cursor_tail_lt b bi hbi hb term hterm
        omega)
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hequal]
    exact ih hnext
  · intro ai bi acc hmore haMore hbMore hnotGreater hnotEqual ih hcursor
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have hless : a[ai]!.1.deg > b[bi]!.1.deg := by omega
    have hamem : a[ai]! ∈ a.toList := by
      rw [getElem!_pos a ai hai]
      exact Array.getElem_mem_toList hai
    have hcoeff : a[ai]!.2 ≠ 0 := ha.2 a[ai]! hamem
    have hnext := hcursor.push a b ai bi (ai + 1) bi acc a[ai]!.1.deg
      a[ai]!.2 hcoeff
      (fun output houtput => hcursor.result_above_a output houtput a[ai]!
        (by rw [drop_eq_getElem_cons a ai hai]; simp))
      (canonical_cursor_tail_lt a ai hai ha)
      (fun term hterm => lt_of_le_of_lt
        (canonical_cursor_degree_le b bi hbi hb term hterm) hless)
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp only [hmore, haMore, hbMore, hnotGreater, hnotEqual, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hdone hcursor
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simpa [hdone] using hcursor.result_canonical

/-- Canonicality of the public generated sparse-add entry. -/
theorem pairVecAdd_canonical (a b : SparsePolyZZ)
    (ha : StrictPolynomialMod.SparsePolyZZCanonical a)
    (hb : StrictPolynomialMod.SparsePolyZZCanonical b) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.pairVecAddLoop a b 0 0 #[]) := by
  apply pairVecAddLoop_canonical a b ha hb 0 0 #[]
  exact ⟨by simp [StrictPolynomialMod.SparsePolyZZCanonical], by simp,
    by simp⟩

/-- The generated sparse merge only appends to its explicit result cursor.
This frame fact follows the actual six-way source recursion and is used to
retain a leading term after the first merge step emits it. -/
theorem pairVecAddLoop_result_prefix (a b : SparsePolyZZ) :
    ∀ aIndex bIndex result,
      ∃ suffix,
        (Generated.StrictHensel.pairVecAddLoop a b aIndex bIndex result).toList =
          result.toList ++ suffix := by
  intro aIndex bIndex result
  refine Generated.StrictHensel.pairVecAddLoop.induct a b
    (motive := fun aIndex bIndex result => ∃ suffix,
      (Generated.StrictHensel.pairVecAddLoop a b aIndex bIndex result).toList =
        result.toList ++ suffix)
    ?_ ?_ ?_ ?_ ?_ ?_ aIndex bIndex result
  · intro ai bi acc hmore haDone ih
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haDone]
    rcases ih with ⟨suffix, hsuffix⟩
    refine ⟨b[bi]! :: suffix, ?_⟩
    rw [hsuffix, Array.toList_push]
    simp
  · intro ai bi acc hmore haMore hbDone ih
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbDone]
    rcases ih with ⟨suffix, hsuffix⟩
    refine ⟨a[ai]! :: suffix, ?_⟩
    rw [hsuffix, Array.toList_push]
    simp
  · intro ai bi acc hmore haMore hbMore hdegree ih
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbMore, hdegree]
    rcases ih with ⟨suffix, hsuffix⟩
    refine ⟨b[bi]! :: suffix, ?_⟩
    rw [hsuffix, Array.toList_push]
    simp
  · intro ai bi acc hmore haMore hbMore hnotGreater hequal ih
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    by_cases hcoefficient : b[bi]!.2 + a[ai]!.2 != 0
    · simp [hmore, haMore, hbMore, hnotGreater, hequal,
        Generated.StrictHensel.pushNonzero, hcoefficient] at ih ⊢
      rcases ih with ⟨suffix, hsuffix⟩
      refine ⟨(a[ai]!.1, b[bi]!.2 + a[ai]!.2) :: suffix, ?_⟩
      simpa using hsuffix
    · simpa [hmore, haMore, hbMore, hnotGreater, hequal,
        Generated.StrictHensel.pushNonzero, hcoefficient] using ih
  · intro ai bi acc hmore haMore hbMore hnotGreater hnotEqual ih
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hnotEqual]
    rcases ih with ⟨suffix, hsuffix⟩
    refine ⟨a[ai]! :: suffix, ?_⟩
    rw [hsuffix, Array.toList_push]
    simp
  · intro ai bi acc hdone
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    exact ⟨[], by simp [hdone]⟩

/-- If every term on the right is below the concrete left head, the actual
merge emits that left head first and the recursive frame keeps it there. -/
theorem pairVecAddLoop_preserves_left_head
    (a b : SparsePolyZZ) (ha : 0 < a.size)
    (hbelow : ∀ term ∈ b.toList, term.1.deg < a[0].1.deg) :
    ∃ suffix,
      (Generated.StrictHensel.pairVecAddLoop a b 0 0 #[]).toList =
        a[0] :: suffix := by
  rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
  have hmore : 0 < a.size ∨ 0 < b.size := Or.inl ha
  have haMore : ¬0 ≥ a.size := by omega
  simp only [hmore, haMore, ↓reduceDIte, ↓reduceIte]
  by_cases hb : 0 < b.size
  · have hbDone : ¬0 ≥ b.size := by omega
    have hbmem : b[0] ∈ b.toList := Array.getElem_mem_toList hb
    have hdegree : ¬b[0]!.1.deg > a[0]!.1.deg := by
      rw [getElem!_pos b 0 hb, getElem!_pos a 0 ha]
      exact Nat.not_lt.mpr (Nat.le_of_lt (hbelow b[0] hbmem))
    have hequal : ¬b[0]!.1.deg = a[0]!.1.deg := by
      rw [getElem!_pos b 0 hb, getElem!_pos a 0 ha]
      exact Nat.ne_of_lt (hbelow b[0] hbmem)
    simp only [hbDone, hdegree, hequal, ↓reduceIte]
    rcases pairVecAddLoop_result_prefix a b 1 0
        ((#[] : SparsePolyZZ).push a[0]!) with
      ⟨suffix, hsuffix⟩
    refine ⟨suffix, ?_⟩
    simpa [getElem!_pos a 0 ha] using hsuffix
  · have hbDone : 0 ≥ b.size := by omega
    simp only [hbDone, ↓reduceIte]
    rcases pairVecAddLoop_result_prefix a b 1 0
        ((#[] : SparsePolyZZ).push a[0]!) with
      ⟨suffix, hsuffix⟩
    refine ⟨suffix, ?_⟩
    simpa [getElem!_pos a 0 ha] using hsuffix

/-- The generated coefficient-scaling range-for changes only coefficients;
every emitted correction term retains its source degree. -/
theorem scaleCoeffs_degrees_below (input : SparsePolyZZ) (scale : Int)
    (bound : Nat) (hbelow : ∀ term ∈ input.toList, term.1.deg < bound) :
    ∀ term ∈ (Generated.StrictHensel.scaleCoeffs input scale).toList,
      term.1.deg < bound := by
  intro term hterm
  rw [Generated.StrictHensel.scaleCoeffs, Array.toList_map,
    List.mem_map] at hterm
  rcases hterm with ⟨source, hsource, rfl⟩
  exact hbelow source hsource

/-- The generated coefficient-reduction/filter loop retains a concrete head
whose floor remainder is itself and whose coefficient is nonzero.  The
statement records the physical output-list head, not merely a modular
polynomial equality. -/
theorem modCoeffOutput_preserves_head
    (input : SparsePolyZZ) (modulus : Int) (head : UMonomial × Int)
    (tail : List (UMonomial × Int))
    (hinput : input.toList = head :: tail)
    (hremainder : ZZ.fdiv_r head.2 head.2 modulus = head.2)
    (hcoefficient : head.2 ≠ 0) :
    ∃ suffix,
      (Generated.StrictHensel.modCoeffOutput input modulus).toList =
        head :: suffix := by
  rw [Generated.StrictHensel.modCoeffOutput, Array.toList_filterMap, hinput]
  simp only [List.filterMap_cons]
  rw [hremainder]
  simp [hcoefficient]

/-- In particular, reduction modulo a positive modulus greater than one
physically retains a leading coefficient equal to one. -/
theorem modCoeffOutput_preserves_one_head
    (input : SparsePolyZZ) (modulus : Nat) (head : UMonomial)
    (tail : List (UMonomial × Int))
    (hinput : input.toList = (head, 1) :: tail)
    (hmodulus : 1 < modulus) :
    ∃ suffix,
      (Generated.StrictHensel.modCoeffOutput input (modulus : Int)).toList =
        (head, 1) :: suffix := by
  apply modCoeffOutput_preserves_head input (modulus : Int) (head, 1) tail
    hinput
  · unfold ZZ.fdiv_r
    change (1 : Int).fmod (modulus : Int) = 1
    exact Int.fmod_eq_of_lt (by omega) (by exact_mod_cast hmodulus)
  · norm_num

/-- The complete concrete `h` correction fragment from `__hensel_step`
preserves a coefficient-one head whenever the actual modular-division
remainder lies below it.  This composes the generated scale, merge, and
coefficient-reduction traversals in source order. -/
theorem henselHCorrection_preserves_one_head
    (h remainder : SparsePolyZZ) (m : Nat) (head : UMonomial)
    (tail : List (UMonomial × Int))
    (hh : h.toList = (head, 1) :: tail)
    (hremainderBelow : ∀ term ∈ remainder.toList,
      term.1.deg < head.deg)
    (hm : 1 < m) :
    ∃ suffix,
      (Generated.StrictHensel.modCoeffOutput
        (Generated.StrictHensel.pairVecAddLoop h
          (Generated.StrictHensel.scaleCoeffs remainder (m : Int))
          0 0 #[]) (m ^ 2 : Int)).toList = (head, 1) :: suffix := by
  have hhsize : 0 < h.size := by
    have hlength := congrArg List.length hh
    have hsizeEq : h.size = tail.length + 1 := by simpa using hlength
    omega
  have hhhead : h[0] = (head, 1) := by
    have hheadList := congrArg List.head? hh
    have hheadOption : h[0]? = some (head, 1) := by
      rw [List.head?_eq_getElem?, Array.getElem?_toList] at hheadList
      exact hheadList
    exact Option.some.inj
      ((Array.getElem?_eq_getElem hhsize).symm.trans hheadOption)
  have hscaledBelow := scaleCoeffs_degrees_below remainder (m : Int)
    head.deg hremainderBelow
  have hadd := pairVecAddLoop_preserves_left_head h
    (Generated.StrictHensel.scaleCoeffs remainder (m : Int)) hhsize (by
      simpa [hhhead] using hscaledBelow)
  rcases hadd with ⟨mergedTail, hmerged⟩
  rw [hhhead] at hmerged
  apply modCoeffOutput_preserves_one_head _ (m ^ 2) head mergedTail hmerged
  nlinarith

set_option maxHeartbeats 0 in
/-- The exact generated `pair_vec_sub` merge preserves the canonical sparse
integer representation under the same emitted-prefix frontier invariant. -/
theorem pairVecSubLoop_canonical (a b : SparsePolyZZ)
    (ha : StrictPolynomialMod.SparsePolyZZCanonical a)
    (hb : StrictPolynomialMod.SparsePolyZZCanonical b) :
    ∀ aIndex bIndex result,
      PairVecAddCursorValid a b aIndex bIndex result →
      StrictPolynomialMod.SparsePolyZZCanonical
        (Generated.StrictHensel.pairVecSubLoop a b aIndex bIndex result) := by
  intro aIndex bIndex result hvalid
  refine Generated.StrictHensel.pairVecSubLoop.induct a b
    (motive := fun aIndex bIndex result =>
      PairVecAddCursorValid a b aIndex bIndex result →
      StrictPolynomialMod.SparsePolyZZCanonical
        (Generated.StrictHensel.pairVecSubLoop a b aIndex bIndex result))
    ?_ ?_ ?_ ?_ ?_ ?_ aIndex bIndex result hvalid
  · intro ai bi acc hmore haDone ih hcursor
    have hbi : bi < b.size := by omega
    have hbmem : b[bi]! ∈ b.toList := by
      rw [getElem!_pos b bi hbi]
      exact Array.getElem_mem_toList hbi
    have hcoeff : -b[bi]!.2 ≠ 0 := neg_ne_zero.mpr (hb.2 b[bi]! hbmem)
    have haDrop : a.toList.drop ai = [] := List.drop_eq_nil_of_le (by
      simpa using haDone)
    have hnext := hcursor.push a b ai bi ai (bi + 1) acc b[bi]!.1.deg
      (-b[bi]!.2) hcoeff
      (fun output houtput => hcursor.result_above_b output houtput b[bi]!
        (by rw [drop_eq_getElem_cons b bi hbi]; simp))
      (by simp [haDrop])
      (canonical_cursor_tail_lt b bi hbi hb)
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp only [hmore, haDone, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hmore haMore hbDone ih hcursor
    have hai : ai < a.size := by omega
    have hamem : a[ai]! ∈ a.toList := by
      rw [getElem!_pos a ai hai]
      exact Array.getElem_mem_toList hai
    have hcoeff : a[ai]!.2 ≠ 0 := ha.2 a[ai]! hamem
    have hbDrop : b.toList.drop bi = [] := List.drop_eq_nil_of_le (by
      simpa using hbDone)
    have hnext := hcursor.push a b ai bi (ai + 1) bi acc a[ai]!.1.deg
      a[ai]!.2 hcoeff
      (fun output houtput => hcursor.result_above_a output houtput a[ai]!
        (by rw [drop_eq_getElem_cons a ai hai]; simp))
      (canonical_cursor_tail_lt a ai hai ha)
      (by simp [hbDrop])
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp only [hmore, haMore, hbDone, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hmore haMore hbMore hdegree ih hcursor
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have hbmem : b[bi]! ∈ b.toList := by
      rw [getElem!_pos b bi hbi]
      exact Array.getElem_mem_toList hbi
    have hcoeff : -b[bi]!.2 ≠ 0 := neg_ne_zero.mpr (hb.2 b[bi]! hbmem)
    have hnext := hcursor.push a b ai bi ai (bi + 1) acc b[bi]!.1.deg
      (-b[bi]!.2) hcoeff
      (fun output houtput => hcursor.result_above_b output houtput b[bi]!
        (by rw [drop_eq_getElem_cons b bi hbi]; simp))
      (fun term hterm => lt_of_le_of_lt
        (canonical_cursor_degree_le a ai hai ha term hterm) hdegree)
      (canonical_cursor_tail_lt b bi hbi hb)
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp only [hmore, haMore, hbMore, hdegree, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hmore haMore hbMore hnotGreater hequal ih hcursor
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have habove : ∀ output ∈ acc.toList,
        a[ai]!.1.deg < output.1.deg := by
      intro output houtput
      exact hcursor.result_above_a output houtput a[ai]!
        (by rw [drop_eq_getElem_cons a ai hai]; simp)
    have hnext := hcursor.pushNonzero a b ai bi (ai + 1) (bi + 1) acc
      a[ai]!.1.deg (a[ai]!.2 - b[bi]!.2) habove
      (canonical_cursor_tail_lt a ai hai ha)
      (fun term hterm => by
        have hlt := canonical_cursor_tail_lt b bi hbi hb term hterm
        omega)
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hequal]
    exact ih hnext
  · intro ai bi acc hmore haMore hbMore hnotGreater hnotEqual ih hcursor
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have hless : a[ai]!.1.deg > b[bi]!.1.deg := by omega
    have hamem : a[ai]! ∈ a.toList := by
      rw [getElem!_pos a ai hai]
      exact Array.getElem_mem_toList hai
    have hcoeff : a[ai]!.2 ≠ 0 := ha.2 a[ai]! hamem
    have hnext := hcursor.push a b ai bi (ai + 1) bi acc a[ai]!.1.deg
      a[ai]!.2 hcoeff
      (fun output houtput => hcursor.result_above_a output houtput a[ai]!
        (by rw [drop_eq_getElem_cons a ai hai]; simp))
      (canonical_cursor_tail_lt a ai hai ha)
      (fun term hterm => lt_of_le_of_lt
        (canonical_cursor_degree_le b bi hbi hb term hterm) hless)
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp only [hmore, haMore, hbMore, hnotGreater, hnotEqual, ↓reduceIte]
    exact ih hnext
  · intro ai bi acc hdone hcursor
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simpa [hdone] using hcursor.result_canonical

/-- Canonicality of the public generated sparse-subtraction entry. -/
theorem pairVecSub_canonical (a b : SparsePolyZZ)
    (ha : StrictPolynomialMod.SparsePolyZZCanonical a)
    (hb : StrictPolynomialMod.SparsePolyZZCanonical b) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.pairVecSubLoop a b 0 0 #[]) := by
  apply pairVecSubLoop_canonical a b ha hb 0 0 #[]
  exact ⟨by simp [StrictPolynomialMod.SparsePolyZZCanonical], by simp,
    by simp⟩

private theorem mulMaxDegree_spec
    (products : List Generated.StrictHensel.MulProduct)
    (hnonempty : products ≠ []) :
    (∃ product ∈ products,
      product.degree = Generated.StrictHensel.mulMaxDegree products) ∧
    ∀ product ∈ products,
      product.degree ≤ Generated.StrictHensel.mulMaxDegree products := by
  have hsome :
      (products.map (fun product => product.degree)).max? =
        some (Generated.StrictHensel.mulMaxDegree products) := by
    unfold Generated.StrictHensel.mulMaxDegree
    cases hmax : (products.map (fun product => product.degree)).max? with
    | none =>
        have hempty := List.max?_eq_none_iff.mp hmax
        simp at hempty
        contradiction
    | some maximum => simp [hmax]
  have hspec := hsome
  rw [List.max?_eq_some_iff] at hspec
  constructor
  · rcases List.mem_map.mp hspec.1 with ⟨product, hproduct, heq⟩
    exact ⟨product, hproduct, heq⟩
  · intro product hproduct
    exact hspec.2 product.degree (List.mem_map.mpr ⟨product, hproduct, rfl⟩)

private structure PairVecMulHeapCursorValid
    (products : List Generated.StrictHensel.MulProduct)
    (result : SparsePolyZZ) : Prop where
  result_canonical : StrictPolynomialMod.SparsePolyZZCanonical result
  result_above_products : ∀ output ∈ result.toList,
    ∀ product ∈ products, product.degree < output.1.deg

private theorem PairVecMulHeapCursorValid.next
    (products : List Generated.StrictHensel.MulProduct)
    (result : SparsePolyZZ) (hnonempty : products ≠ [])
    (hvalid : PairVecMulHeapCursorValid products result) :
    let degree := Generated.StrictHensel.mulMaxDegree products
    let coefficient := Generated.StrictHensel.mulDegreeCoefficient degree products
    let remaining := products.filter fun product => product.degree != degree
    PairVecMulHeapCursorValid remaining
      (Generated.StrictHensel.pushNonzero result degree coefficient) := by
  dsimp only
  have hspec := mulMaxDegree_spec products hnonempty
  let degree := Generated.StrictHensel.mulMaxDegree products
  let coefficient := Generated.StrictHensel.mulDegreeCoefficient degree products
  let remaining := products.filter fun product => product.degree != degree
  have habove : ∀ output ∈ result.toList, degree < output.1.deg := by
    intro output houtput
    rcases hspec.1 with ⟨product, hproduct, hdegree⟩
    have hlt := hvalid.result_above_products output houtput product hproduct
    simpa [degree, hdegree] using hlt
  refine ⟨pushNonzero_canonical result degree coefficient
      hvalid.result_canonical habove, ?_⟩
  intro output houtput product hproduct
  have hproductSource : product ∈ products := by
    exact List.mem_of_mem_filter hproduct
  have hproductNe : product.degree ≠ degree := by
    simpa [degree] using List.of_mem_filter hproduct
  have hproductLt : product.degree < degree :=
    lt_of_le_of_ne (hspec.2 product hproductSource) hproductNe
  by_cases hzero : coefficient = 0
  · have hzero' : Generated.StrictHensel.mulDegreeCoefficient
        (Generated.StrictHensel.mulMaxDegree products) products = 0 := by
      simpa [coefficient, degree] using hzero
    have hold : output ∈ result.toList := by
      simpa [Generated.StrictHensel.pushNonzero, hzero'] using houtput
    exact lt_trans hproductLt (habove output hold)
  · have hzero' : Generated.StrictHensel.mulDegreeCoefficient
        (Generated.StrictHensel.mulMaxDegree products) products ≠ 0 := by
      simpa [coefficient, degree] using hzero
    rw [Generated.StrictHensel.pushNonzero, if_pos (by simpa using hzero'),
      Array.toList_push] at houtput
    simp only [List.mem_append, List.mem_singleton] at houtput
    rcases houtput with hold | rfl
    · exact lt_trans hproductLt (habove output hold)
    · exact hproductLt

set_option maxHeartbeats 0 in
/-- The exact generated multiplication heap traversal emits a canonical
sparse integer polynomial from every valid emitted-prefix frontier. -/
theorem pairVecMulHeapLoop_canonical :
    ∀ (products : List Generated.StrictHensel.MulProduct)
      (result : SparsePolyZZ),
      PairVecMulHeapCursorValid products result →
      StrictPolynomialMod.SparsePolyZZCanonical
        (Generated.StrictHensel.pairVecMulHeapLoop products result) := by
  intro products result hvalid
  refine Generated.StrictHensel.pairVecMulHeapLoop.induct
    (motive := fun products result =>
      PairVecMulHeapCursorValid products result →
      StrictPolynomialMod.SparsePolyZZCanonical
        (Generated.StrictHensel.pairVecMulHeapLoop products result))
    ?_ ?_ products result hvalid
  · intro result hcursor
    rw [Generated.StrictHensel.pairVecMulHeapLoop.eq_1]
    simpa using hcursor.result_canonical
  · intro products result hempty
    dsimp only
    intro ih hcursor
    have hnext := PairVecMulHeapCursorValid.next products result hempty hcursor
    have hremaining :
        (products.attach.filter fun product =>
          product.1.degree != Generated.StrictHensel.mulMaxDegree products).unattach =
          products.filter fun product =>
            product.degree != Generated.StrictHensel.mulMaxDegree products := by
      simp
    rw [hremaining] at ih
    rw [Generated.StrictHensel.pairVecMulHeapLoop.eq_1]
    simp only [hempty, ↓reduceDIte]
    exact ih hnext

/-- Canonicality of the public generated sparse-multiplication entry. -/
theorem pairVecMul_canonical (a b : SparsePolyZZ) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts a b) #[]) := by
  apply pairVecMulHeapLoop_canonical _ #[]
  exact ⟨by simp [StrictPolynomialMod.SparsePolyZZCanonical], by simp⟩

theorem strictHenselRawOps_mul_canonical_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (a b output : SparsePolyZZ)
    (hrun : (strictHenselRawOps termination).mul a b = .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  simp only [strictHenselRawOps,
    Generated.StrictHensel.__upoly_mul_raw_ir] at hrun
  injection hrun with houtput
  subst output
  exact pairVecMul_canonical a b

theorem strictHenselRawOps_add_canonical_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (a b output : SparsePolyZZ)
    (ha : StrictPolynomialMod.SparsePolyZZCanonical a)
    (hb : StrictPolynomialMod.SparsePolyZZCanonical b)
    (hrun : (strictHenselRawOps termination).add a b = .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  simp only [strictHenselRawOps,
    Generated.StrictHensel.__upoly_add_raw_ir] at hrun
  injection hrun with houtput
  subst output
  exact pairVecAdd_canonical a b ha hb

theorem strictHenselRawOps_sub_canonical_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (a b output : SparsePolyZZ)
    (ha : StrictPolynomialMod.SparsePolyZZCanonical a)
    (hb : StrictPolynomialMod.SparsePolyZZCanonical b)
    (hrun : (strictHenselRawOps termination).sub a b = .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  simp only [strictHenselRawOps,
    Generated.StrictHensel.__upoly_sub_raw_ir] at hrun
  injection hrun with houtput
  subst output
  exact pairVecSub_canonical a b ha hb

theorem modCoeffOutput_canonical_of_run
    (f output : SparsePolyZZ) (m : Int)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hrun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir f m =
      .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  simp only [Generated.StrictHensel.__upoly_mod_coeff_raw_ir] at hrun
  injection hrun with houtput
  subst output
  exact modCoeffOutput_canonical f m hf

private structure DivmodMergeCursorValid (r g : SparsePolyZZ)
    (degreeShift rIndex gIndex : Nat) (result : SparsePolyZZ) : Prop where
  result_canonical : StrictPolynomialMod.SparsePolyZZCanonical result
  result_above_r : ∀ output ∈ result.toList,
    ∀ term ∈ r.toList.drop rIndex, term.1.deg < output.1.deg
  result_above_g : ∀ output ∈ result.toList,
    ∀ term ∈ g.toList.drop gIndex,
      term.1.deg + degreeShift < output.1.deg

private theorem DivmodMergeCursorValid.push
    (r g : SparsePolyZZ) (degreeShift rIndex gIndex : Nat)
    (result : SparsePolyZZ) (degree : Nat) (coefficient : Int)
    (hresult : StrictPolynomialMod.SparsePolyZZCanonical result)
    (habove : ∀ output ∈ result.toList, degree < output.1.deg)
    (hrBelow : ∀ term ∈ r.toList.drop rIndex, term.1.deg < degree)
    (hgBelow : ∀ term ∈ g.toList.drop gIndex,
      term.1.deg + degreeShift < degree) :
    DivmodMergeCursorValid r g degreeShift rIndex gIndex
      (Generated.StrictHensel.pushNonzero result degree coefficient) := by
  refine ⟨pushNonzero_canonical result degree coefficient hresult habove,
    ?_, ?_⟩
  · intro output houtput term hterm
    by_cases hzero : coefficient = 0
    · subst coefficient
      have hout : output ∈ result.toList := by
        simpa [Generated.StrictHensel.pushNonzero] using houtput
      exact lt_trans (hrBelow term hterm) (habove output hout)
    · unfold Generated.StrictHensel.pushNonzero at houtput
      simp only [bne_iff_ne] at houtput
      rw [if_pos hzero, Array.toList_push] at houtput
      simp only [List.mem_append, List.mem_singleton] at houtput
      rcases houtput with hold | rfl
      · exact lt_trans (hrBelow term hterm) (habove output hold)
      · exact hrBelow term hterm
  · intro output houtput term hterm
    by_cases hzero : coefficient = 0
    · subst coefficient
      have hout : output ∈ result.toList := by
        simpa [Generated.StrictHensel.pushNonzero] using houtput
      exact lt_trans (hgBelow term hterm) (habove output hout)
    · unfold Generated.StrictHensel.pushNonzero at houtput
      simp only [bne_iff_ne] at houtput
      rw [if_pos hzero, Array.toList_push] at houtput
      simp only [List.mem_append, List.mem_singleton] at houtput
      rcases houtput with hold | rfl
      · exact lt_trans (hgBelow term hterm) (habove output hold)
      · exact hgBelow term hterm

private theorem DivmodMergeCursorValid.takeR
    (r g : SparsePolyZZ) (degreeShift rIndex gIndex : Nat)
    (result : SparsePolyZZ) (coefficient : Int)
    (hr : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hrIndex : rIndex < r.size)
    (hcursor : DivmodMergeCursorValid r g degreeShift rIndex gIndex result)
    (hgBelow : ∀ term ∈ g.toList.drop gIndex,
      term.1.deg + degreeShift < r[rIndex]!.1.deg) :
    DivmodMergeCursorValid r g degreeShift (rIndex + 1) gIndex
      (Generated.StrictHensel.pushNonzero result r[rIndex]!.1.deg
        coefficient) := by
  apply DivmodMergeCursorValid.push r g degreeShift (rIndex + 1) gIndex
    result r[rIndex]!.1.deg coefficient hcursor.result_canonical
  · intro output houtput
    apply hcursor.result_above_r output houtput r[rIndex]!
    rw [sparseDrop_eq_getElem_cons r rIndex hrIndex]
    simp
  · exact canonical_cursor_tail_lt r rIndex hrIndex hr
  · exact hgBelow

private theorem DivmodMergeCursorValid.takeG
    (r g : SparsePolyZZ) (degreeShift rIndex gIndex : Nat)
    (result : SparsePolyZZ) (coefficient : Int)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hgIndex : gIndex < g.size)
    (hcursor : DivmodMergeCursorValid r g degreeShift rIndex gIndex result)
    (hrBelow : ∀ term ∈ r.toList.drop rIndex,
      term.1.deg < g[gIndex]!.1.deg + degreeShift) :
    DivmodMergeCursorValid r g degreeShift rIndex (gIndex + 1)
      (Generated.StrictHensel.pushNonzero result
        (g[gIndex]!.1.deg + degreeShift) coefficient) := by
  apply DivmodMergeCursorValid.push r g degreeShift rIndex (gIndex + 1)
    result (g[gIndex]!.1.deg + degreeShift) coefficient
    hcursor.result_canonical
  · intro output houtput
    apply hcursor.result_above_g output houtput g[gIndex]!
    rw [sparseDrop_eq_getElem_cons g gIndex hgIndex]
    simp
  · exact hrBelow
  · intro term hterm
    have hlt := canonical_cursor_tail_lt g gIndex hgIndex hg term hterm
    omega

private theorem DivmodMergeCursorValid.takeBoth
    (r g : SparsePolyZZ) (degreeShift rIndex gIndex : Nat)
    (result : SparsePolyZZ) (coefficient : Int)
    (hr : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hrIndex : rIndex < r.size) (hgIndex : gIndex < g.size)
    (hcursor : DivmodMergeCursorValid r g degreeShift rIndex gIndex result)
    (hdegree : r[rIndex]!.1.deg = g[gIndex]!.1.deg + degreeShift) :
    DivmodMergeCursorValid r g degreeShift (rIndex + 1) (gIndex + 1)
      (Generated.StrictHensel.pushNonzero result r[rIndex]!.1.deg
        coefficient) := by
  apply DivmodMergeCursorValid.push r g degreeShift (rIndex + 1)
    (gIndex + 1) result r[rIndex]!.1.deg coefficient
    hcursor.result_canonical
  · intro output houtput
    apply hcursor.result_above_r output houtput r[rIndex]!
    rw [sparseDrop_eq_getElem_cons r rIndex hrIndex]
    simp
  · exact canonical_cursor_tail_lt r rIndex hrIndex hr
  · intro term hterm
    have hlt := canonical_cursor_tail_lt g gIndex hgIndex hg term hterm
    omega

theorem divmodMergeLoop_canonical (r g : SparsePolyZZ)
    (coefficient m : Int) (degreeShift rIndex gIndex : Nat)
    (result : SparsePolyZZ)
    (hr : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hcursor : DivmodMergeCursorValid r g degreeShift rIndex gIndex result) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.divmodMergeLoop r g coefficient m degreeShift
        rIndex gIndex result) := by
  refine Generated.StrictHensel.divmodMergeLoop.induct r g coefficient m
    degreeShift
    (motive := fun ri gi acc =>
      DivmodMergeCursorValid r g degreeShift ri gi acc →
        StrictPolynomialMod.SparsePolyZZCanonical
          (Generated.StrictHensel.divmodMergeLoop r g coefficient m
            degreeShift ri gi acc)) ?_ ?_ ?_ ?_ ?_ ?_ rIndex gIndex result
    hcursor
  · intro ri gi acc hmore hgDone
    dsimp only
    intro ih hvalid
    have hri : ri < r.size := by omega
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp only [hmore, hgDone, ↓reduceDIte, ↓reduceIte]
    apply ih
    apply DivmodMergeCursorValid.takeR r g degreeShift ri gi acc _ hr hri
      hvalid
    intro term hterm
    have hgNil : g.toList.drop gi = [] :=
      List.drop_eq_nil_iff.mpr (by simpa using hgDone)
    rw [hgNil] at hterm
    contradiction
  · intro ri gi acc hmore hgMore hrDone
    dsimp only
    intro ih hvalid
    have hgi : gi < g.size := by omega
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp only [hmore, hgMore, hrDone, ↓reduceDIte, ↓reduceIte]
    apply ih
    apply DivmodMergeCursorValid.takeG r g degreeShift ri gi acc _ hg hgi
      hvalid
    intro term hterm
    have hrNil : r.toList.drop ri = [] :=
      List.drop_eq_nil_iff.mpr (by simpa using hrDone)
    rw [hrNil] at hterm
    contradiction
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hdegree ih hvalid
    have hri : ri < r.size := by omega
    have hgi : gi < g.size := by omega
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp only [hmore, hgMore, hrMore, hdegree, ↓reduceDIte, ↓reduceIte]
    apply ih
    apply DivmodMergeCursorValid.takeR r g degreeShift ri gi acc _ hr hri
      hvalid
    intro term hterm
    have hle := canonical_cursor_degree_le g gi hgi hg term hterm
    omega
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hless ih hvalid
    have hri : ri < r.size := by omega
    have hgi : gi < g.size := by omega
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp only [hmore, hgMore, hrMore, hnotGreater, hless, ↓reduceDIte,
      ↓reduceIte]
    apply ih
    apply DivmodMergeCursorValid.takeG r g degreeShift ri gi acc _ hg hgi
      hvalid
    intro term hterm
    have hle := canonical_cursor_degree_le r ri hri hr term hterm
    omega
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hnotLess ih hvalid
    have hri : ri < r.size := by omega
    have hgi : gi < g.size := by omega
    have hequal : r[ri]!.1.deg = g[gi]!.1.deg + degreeShift := by omega
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp only [hmore, hgMore, hrMore, hnotGreater, hnotLess, ↓reduceDIte,
      ↓reduceIte]
    apply ih
    exact DivmodMergeCursorValid.takeBoth r g degreeShift ri gi acc _ hr hg
      hri hgi hvalid hequal
  · intro ri gi acc hdone hvalid
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simpa [hdone] using hvalid.result_canonical

theorem divmodRemainder_canonical (r g : SparsePolyZZ)
    (coefficient m : Int) (degreeShift : Nat)
    (hr : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.divmodRemainder r g coefficient m degreeShift) := by
  unfold Generated.StrictHensel.divmodRemainder
  apply divmodMergeLoop_canonical r g coefficient m degreeShift 1 1 #[] hr hg
  exact ⟨by simp [StrictPolynomialMod.SparsePolyZZCanonical], by simp,
    by simp⟩

/-- Every exact finite trace consumed by the generated modular long-division
loop returns a canonical remainder when its current remainder and fixed
divisor are canonical. -/
theorem Generated.StrictHensel.DivmodTrace.remainder_canonical
    {g : SparsePolyZZ} {inverse m : ZZ} {r q : SparsePolyZZ}
    (trace : Generated.StrictHensel.DivmodTrace g inverse m r q)
    (hr : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.divmodLoop g inverse m trace).2 := by
  induction trace with
  | done r q inactive =>
      simpa [Generated.StrictHensel.divmodLoop] using hr
  | vanished r q active zero next ih =>
      have hparts : r.isEmpty = false ∧ get_deg r ≥ get_deg g := by
        simpa using active
      have hrsize : 0 < r.size := by
        have hne : r ≠ #[] := by
          intro hempty
          subst r
          simp at hparts
        have hsizeNe : r.size ≠ 0 := by
          intro hzero
          exact hne (Array.eq_empty_of_size_eq_zero hzero)
        omega
      rw [Generated.StrictHensel.divmodLoop]
      exact ih (eraseLeading_canonical r hrsize hr)
  | subtract r q active nonzero next ih =>
      rw [Generated.StrictHensel.divmodLoop]
      exact ih (divmodRemainder_canonical r g
        (Generated.StrictHensel.divmodCoefficient r inverse m) m
        (get_deg r - get_deg g).toNatClampNeg hr hg)

private def DegreesBelow (bound : Nat) (f : SparsePolyZZ) : Prop :=
  ∀ term ∈ f.toList, term.1.deg < bound

private theorem degreesBelow_pushNonzero (bound : Nat) (f : SparsePolyZZ)
    (degree : Nat) (coefficient : Int) (hf : DegreesBelow bound f)
    (hdegree : degree < bound) :
    DegreesBelow bound
      (Generated.StrictHensel.pushNonzero f degree coefficient) := by
  intro term hterm
  by_cases hzero : coefficient = 0
  · subst coefficient
    apply hf term
    simpa [Generated.StrictHensel.pushNonzero] using hterm
  · unfold Generated.StrictHensel.pushNonzero at hterm
    simp only [bne_iff_ne] at hterm
    rw [if_pos hzero, Array.toList_push] at hterm
    simp only [List.mem_append, List.mem_singleton] at hterm
    rcases hterm with hold | rfl
    · exact hf term hold
    · exact hdegree

private theorem mem_drop_of_mem_drop_succ {alpha : Type*} (input : List alpha)
    (index : Nat) {value : alpha} (hvalue : value ∈ input.drop (index + 1)) :
    value ∈ input.drop index := by
  have heq : input.drop (index + 1) = (input.drop index).drop 1 := by
    rw [List.drop_drop]
  rw [heq] at hvalue
  exact List.mem_of_mem_drop hvalue

theorem divmodMergeLoop_degreesBelow (bound : Nat) (r g : SparsePolyZZ)
    (coefficient m : Int) (degreeShift : Nat) :
    ∀ (rIndex gIndex : Nat) (result : SparsePolyZZ),
      DegreesBelow bound result →
      (∀ term ∈ r.toList.drop rIndex, term.1.deg < bound) →
      (∀ term ∈ g.toList.drop gIndex,
        term.1.deg + degreeShift < bound) →
      DegreesBelow bound
        (Generated.StrictHensel.divmodMergeLoop r g coefficient m degreeShift
          rIndex gIndex result) := by
  intro rIndex gIndex result
  refine Generated.StrictHensel.divmodMergeLoop.induct r g coefficient m
    degreeShift
    (motive := fun ri gi acc => DegreesBelow bound acc →
      (∀ term ∈ r.toList.drop ri, term.1.deg < bound) →
      (∀ term ∈ g.toList.drop gi,
        term.1.deg + degreeShift < bound) →
      DegreesBelow bound
        (Generated.StrictHensel.divmodMergeLoop r g coefficient m degreeShift
          ri gi acc)) ?_ ?_ ?_ ?_ ?_ ?_ rIndex gIndex result
  · intro ri gi acc hmore hgDone
    dsimp only
    intro ih hacc hr hg
    have hri : ri < r.size := by omega
    have hcurrent : r[ri]!.1.deg < bound := by
      apply hr r[ri]!
      rw [sparseDrop_eq_getElem_cons r ri hri]
      simp
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgDone]
    exact ih (degreesBelow_pushNonzero bound acc r[ri]!.1.deg _ hacc hcurrent)
      (fun term hterm => hr term
        (mem_drop_of_mem_drop_succ r.toList ri hterm)) hg
  · intro ri gi acc hmore hgMore hrDone
    dsimp only
    intro ih hacc hr hg
    have hgi : gi < g.size := by omega
    have hcurrent : g[gi]!.1.deg + degreeShift < bound := by
      apply hg g[gi]!
      rw [sparseDrop_eq_getElem_cons g gi hgi]
      simp
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrDone]
    simp only [Int.emod_emod] at ih
    exact ih (degreesBelow_pushNonzero bound acc
      (g[gi]!.1.deg + degreeShift) _ hacc hcurrent) hr
      (fun term hterm => hg term
        (mem_drop_of_mem_drop_succ g.toList gi hterm))
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hdegree ih hacc hr hg
    have hri : ri < r.size := by omega
    have hcurrent : r[ri]!.1.deg < bound := by
      apply hr r[ri]!
      rw [sparseDrop_eq_getElem_cons r ri hri]
      simp
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hdegree]
    exact ih (degreesBelow_pushNonzero bound acc r[ri]!.1.deg _ hacc hcurrent)
      (fun term hterm => hr term
        (mem_drop_of_mem_drop_succ r.toList ri hterm)) hg
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hless ih hacc hr hg
    have hgi : gi < g.size := by omega
    have hcurrent : g[gi]!.1.deg + degreeShift < bound := by
      apply hg g[gi]!
      rw [sparseDrop_eq_getElem_cons g gi hgi]
      simp
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hnotGreater, hless]
    simp only [Int.emod_emod] at ih
    exact ih (degreesBelow_pushNonzero bound acc
      (g[gi]!.1.deg + degreeShift) _ hacc hcurrent) hr
      (fun term hterm => hg term
        (mem_drop_of_mem_drop_succ g.toList gi hterm))
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hnotLess ih hacc hr hg
    have hri : ri < r.size := by omega
    have hcurrent : r[ri]!.1.deg < bound := by
      apply hr r[ri]!
      rw [sparseDrop_eq_getElem_cons r ri hri]
      simp
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hnotGreater, hnotLess]
    exact ih (degreesBelow_pushNonzero bound acc r[ri]!.1.deg _ hacc hcurrent)
      (fun term hterm => hr term
        (mem_drop_of_mem_drop_succ r.toList ri hterm))
      (fun term hterm => hg term
        (mem_drop_of_mem_drop_succ g.toList gi hterm))
  · intro ri gi acc hdone hacc hr hg
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simpa [hdone] using hacc

theorem divmodRemainder_degreesBelow_head (r g : SparsePolyZZ)
    (coefficient m : Int) (hrsize : 0 < r.size) (hgsize : 0 < g.size)
    (hrCanonical : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hgCanonical : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hrDegree : r[0]!.1.deg < 2 ^ 63)
    (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hactive : get_deg r ≥ get_deg g) :
    DegreesBelow r[0]!.1.deg
      (Generated.StrictHensel.divmodRemainder r g coefficient m
        (get_deg r - get_deg g).toNatClampNeg) := by
  have hshift := get_deg_sub_toNatClampNeg_eq_shift r g hrsize hgsize
    hrDegree hgDegree hactive
  unfold Generated.StrictHensel.divmodRemainder
  apply divmodMergeLoop_degreesBelow r[0]!.1.deg r g coefficient m
    (get_deg r - get_deg g).toNatClampNeg 1 1 #[]
  · simp [DegreesBelow]
  · exact canonical_cursor_tail_lt r 0 hrsize hrCanonical
  · intro term hterm
    have hlt := canonical_cursor_tail_lt g 0 hgsize hgCanonical term hterm
    omega

theorem divmodRemainder_head_lt (r g output : SparsePolyZZ)
    (coefficient m : Int) (hrsize : 0 < r.size) (hgsize : 0 < g.size)
    (houtputSize : 0 < output.size)
    (hrCanonical : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hgCanonical : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hrDegree : r[0]!.1.deg < 2 ^ 63)
    (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hactive : get_deg r ≥ get_deg g)
    (houtput : Generated.StrictHensel.divmodRemainder r g coefficient m
      (get_deg r - get_deg g).toNatClampNeg = output) :
    output[0]!.1.deg < r[0]!.1.deg := by
  have hbelow := divmodRemainder_degreesBelow_head r g coefficient m hrsize
    hgsize hrCanonical hgCanonical hrDegree hgDegree hactive
  rw [houtput] at hbelow
  apply hbelow output[0]!
  rw [getElem!_pos output 0 houtputSize]
  exact Array.getElem_mem_toList houtputSize

private def divmodDegreeRank (r : SparsePolyZZ) : Nat :=
  if r.isEmpty then 0 else r[0]!.1.deg + 1

theorem eraseLeading_divmodDegreeRank_lt (r : SparsePolyZZ)
    (hrsize : 0 < r.size)
    (hrCanonical : StrictPolynomialMod.SparsePolyZZCanonical r) :
    divmodDegreeRank (r.eraseIdxIfInBounds 0) < divmodDegreeRank r := by
  have hrNonempty : r.isEmpty = false := by
    apply Bool.eq_false_of_not_eq_true
    rw [Array.isEmpty_iff]
    intro hempty
    subst r
    simp at hrsize
  unfold divmodDegreeRank
  rw [hrNonempty]
  simp only [Bool.false_eq_true, if_false]
  by_cases hnext : (r.eraseIdxIfInBounds 0).isEmpty
  · rw [hnext]
    simp
  · rw [show (r.eraseIdxIfInBounds 0).isEmpty = false by
        exact Bool.eq_false_of_not_eq_true hnext]
    simp only [Bool.false_eq_true, if_false]
    have hnextSize : 0 < (r.eraseIdxIfInBounds 0).size := by
      have hne : r.eraseIdxIfInBounds 0 ≠ #[] := by
        simpa [Array.isEmpty_iff] using hnext
      have hsizeNe : (r.eraseIdxIfInBounds 0).size ≠ 0 := by
        intro hzero
        exact hne (Array.eq_empty_of_size_eq_zero hzero)
      omega
    have hcanonical := eraseLeading_canonical r hrsize hrCanonical
    have hhead := canonical_tail_degree_lt_head r hrsize hrCanonical
      (r.eraseIdxIfInBounds 0)[0]! (by
        have hmem : (r.eraseIdxIfInBounds 0)[0]! ∈
            (r.eraseIdxIfInBounds 0).toList := by
          rw [getElem!_pos _ 0 hnextSize]
          exact Array.getElem_mem_toList hnextSize
        rw [Array.toList_eraseIdxIfInBounds, List.eraseIdx_zero] at hmem
        simpa only [List.drop_one] using hmem)
    exact Nat.add_lt_add_right hhead 1

theorem divmodRemainder_divmodDegreeRank_lt (r g : SparsePolyZZ)
    (coefficient m : Int) (hrsize : 0 < r.size) (hgsize : 0 < g.size)
    (hrCanonical : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hgCanonical : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hrDegree : r[0]!.1.deg < 2 ^ 63)
    (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hactive : get_deg r ≥ get_deg g) :
    divmodDegreeRank
        (Generated.StrictHensel.divmodRemainder r g coefficient m
          (get_deg r - get_deg g).toNatClampNeg) <
      divmodDegreeRank r := by
  let output := Generated.StrictHensel.divmodRemainder r g coefficient m
    (get_deg r - get_deg g).toNatClampNeg
  have hrNonempty : r.isEmpty = false := by
    apply Bool.eq_false_of_not_eq_true
    rw [Array.isEmpty_iff]
    intro hempty
    subst r
    simp at hrsize
  unfold divmodDegreeRank
  rw [hrNonempty]
  simp only [Bool.false_eq_true, if_false]
  by_cases houtputEmpty : output.isEmpty
  · rw [show output.isEmpty = true from houtputEmpty]
    simp
  · have houtputFalse : output.isEmpty = false :=
      Bool.eq_false_of_not_eq_true houtputEmpty
    rw [houtputFalse]
    simp only [Bool.false_eq_true, if_false]
    have houtputSize : 0 < output.size := by
      have hne : output ≠ #[] := by
        simpa [Array.isEmpty_iff] using houtputEmpty
      have hsizeNe : output.size ≠ 0 := by
        intro hzero
        exact hne (Array.eq_empty_of_size_eq_zero hzero)
      omega
    have hhead := divmodRemainder_head_lt r g output coefficient m hrsize
      hgsize houtputSize hrCanonical hgCanonical hrDegree hgDegree hactive rfl
    exact Nat.add_lt_add_right hhead 1

theorem canonical_headDominates (f : SparsePolyZZ) (hfsize : 0 < f.size)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f) : HeadDominates f := by
  intro term hterm
  have hle := canonical_cursor_degree_le f 0 hfsize hf term (by
    simpa using hterm)
  exact hle

/-- Concrete well-founded trace of the C++ modular long-division while-loop.
The proof-erased recursion follows the source active/coefficient branches and
uses the current remainder's leading degree as its well-founded rank. -/
def concreteDivmodTrace (g : SparsePolyZZ) (inverse m : Int)
    (hgsize : 0 < g.size)
    (hgCanonical : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (r q : SparsePolyZZ)
    (hrCanonical : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hrBound : DegreesBound r) :
    Generated.StrictHensel.DivmodTrace g inverse m r q :=
  if hactive : !r.isEmpty && get_deg r ≥ get_deg g then
    have hparts : r.isEmpty = false ∧ get_deg r ≥ get_deg g := by
      simpa using hactive
    have hrsize : 0 < r.size := by
      have hne : r ≠ #[] := by
        intro hempty
        subst r
        simp at hparts
      have hsizeNe : r.size ≠ 0 := by
        intro hzero
        exact hne (Array.eq_empty_of_size_eq_zero hzero)
      omega
    if hzero : Generated.StrictHensel.divmodCoefficient r inverse m = 0 then
      .vanished r q hactive hzero
        (concreteDivmodTrace g inverse m hgsize hgCanonical hgDegree
          (r.eraseIdxIfInBounds 0) q
          (eraseLeading_canonical r hrsize hrCanonical)
          (degreesBound_eraseIdxIfInBounds r 0 hrBound))
    else
      let next := Generated.StrictHensel.divmodRemainder r g
        (Generated.StrictHensel.divmodCoefficient r inverse m) m
        (get_deg r - get_deg g).toNatClampNeg
      .subtract r q hactive hzero
        (concreteDivmodTrace g inverse m hgsize hgCanonical hgDegree next
          (q.push (UMonomial.mk (get_deg r - get_deg g).toNatClampNeg,
            Generated.StrictHensel.divmodCoefficient r inverse m))
          (divmodRemainder_canonical r g
            (Generated.StrictHensel.divmodCoefficient r inverse m) m
            (get_deg r - get_deg g).toNatClampNeg hrCanonical hgCanonical)
          (divmodRemainder_degreesBound r g
            (Generated.StrictHensel.divmodCoefficient r inverse m) m
            hrsize
            hgsize hrBound (canonical_headDominates g hgsize hgCanonical)
            (by
              exact hrBound r[0]! (by
                rw [getElem!_pos r 0 hrsize]
                exact Array.getElem_mem_toList hrsize))
            hgDegree hparts.2))
  else .done r q hactive
termination_by divmodDegreeRank r
decreasing_by
  · apply eraseLeading_divmodDegreeRank_lt r
    · exact hrsize
    · exact hrCanonical
  · apply divmodRemainder_divmodDegreeRank_lt r g
      (Generated.StrictHensel.divmodCoefficient r inverse m) m
    · exact hrsize
    · exact hgsize
    · exact hrCanonical
    · exact hgCanonical
    · exact hrBound r[0]! (by
        rw [getElem!_pos r 0 hrsize]
        exact Array.getElem_mem_toList hrsize)
    · exact hgDegree
    · exact hparts.2

private def sparseCanonicalCheck : List (UMonomial × Int) → Bool
  | [] => true
  | term :: [] => term.2 != 0
  | term :: next :: rest =>
      term.2 != 0 && term.1.deg > next.1.deg &&
        sparseCanonicalCheck (next :: rest)

private theorem sparseCanonicalCheck_eq_true (xs : List (UMonomial × Int)) :
    sparseCanonicalCheck xs = true ↔
      List.IsChain (fun a b : UMonomial × Int => a.1.deg > b.1.deg) xs ∧
      ∀ term ∈ xs, term.2 ≠ 0 := by
  induction xs with
  | nil => simp [sparseCanonicalCheck]
  | cons term xs ih =>
      cases xs with
      | nil => simp [sparseCanonicalCheck]
      | cons next rest =>
          simp [sparseCanonicalCheck, ih, and_assoc, and_left_comm, and_comm]

private def degreesBoundCheck (xs : List (UMonomial × Int)) : Bool :=
  xs.all fun term => decide (term.1.deg < 2 ^ 63)

private theorem degreesBoundCheck_eq_true (xs : List (UMonomial × Int)) :
    degreesBoundCheck xs = true ↔
      ∀ term ∈ xs, term.1.deg < 2 ^ 63 := by
  simp [degreesBoundCheck]

/-- Executable check of exactly the source representation domain on which
modular long division is called.  It checks a nonempty canonical divisor and
a canonical reduced dividend, with all exponents representable by the signed
C++ degree type.  It contains no quotient, remainder, trace, or semantic
polynomial result. -/
def concreteDivmodInputValidCheck (g r : SparsePolyZZ) : Bool :=
  decide (0 < g.size) &&
  sparseCanonicalCheck g.toList && degreesBoundCheck g.toList &&
  sparseCanonicalCheck r.toList && degreesBoundCheck r.toList

def ConcreteDivmodInputValid (g r : SparsePolyZZ) : Prop :=
  concreteDivmodInputValidCheck g r = true

theorem concreteDivmodInputValid_iff (g r : SparsePolyZZ) :
    ConcreteDivmodInputValid g r ↔
      0 < g.size ∧
      StrictPolynomialMod.SparsePolyZZCanonical g ∧
      DegreesBound g ∧
      StrictPolynomialMod.SparsePolyZZCanonical r ∧
      DegreesBound r := by
  simp [ConcreteDivmodInputValid, concreteDivmodInputValidCheck,
    sparseCanonicalCheck_eq_true, degreesBoundCheck_eq_true,
    StrictPolynomialMod.SparsePolyZZCanonical, DegreesBound]
  tauto

/-- The unique concrete termination provider used by the strict Hensel and
recombination executions.  Its trace is constructed by well-founded recursion
over the actual generated C++ remainder update. -/
def concreteDivmodTermination : Generated.StrictHensel.DivmodTermination where
  inputValid := ConcreteDivmodInputValid
  inputValidDecidable := by
    intro g r
    unfold ConcreteDivmodInputValid
    infer_instance
  trace := by
    intro f g reduced m _ hvalid
    rw [concreteDivmodInputValid_iff] at hvalid
    rcases hvalid with
      ⟨hgsize, hgCanonical, hgBound, hrCanonical, hrBound⟩
    have hgDegree : g[0]!.1.deg < 2 ^ 63 := by
      apply hgBound g[0]!
      rw [getElem!_pos g 0 hgsize]
      exact Array.getElem_mem_toList hgsize
    exact concreteDivmodTrace g (ZZ.invert 0 g[0]!.2 m).2 m
      hgsize hgCanonical hgDegree reduced #[] hrCanonical hrBound

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

/-- Every transition of the exact generated division trace preserves the
machine exponent bound, so the physical remainder returned by `divmodLoop`
inherits it from the concrete reduced input. -/
theorem divmodLoop_remainder_degreesBound (m : Nat) (g : SparsePolyZZ)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g) :
    ∀ {r q : SparsePolyZZ}
      (trace : Generated.StrictHensel.DivmodTrace g
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) r q),
      DegreesBound r →
      DegreesBound (Generated.StrictHensel.divmodLoop g
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) trace).2 := by
  intro r q trace hrBound
  induction trace with
  | done r q inactive => exact hrBound
  | vanished r q active zero next ih =>
      exact ih (degreesBound_eraseIdxIfInBounds r 0 hrBound)
  | subtract r q active nonzero next ih =>
      have hactive : r.isEmpty = false ∧ get_deg r ≥ get_deg g := by
        simpa using active
      have hr : 0 < r.size := by
        have hrne : r ≠ #[] := by
          simpa [Array.isEmpty_iff] using hactive.1
        have hrsize : r.size ≠ 0 := by
          intro hzero
          exact hrne (Array.eq_empty_of_size_eq_zero hzero)
        omega
      have hrDegree : r[0]!.1.deg < 2 ^ 63 := by
        apply hrBound r[0]!
        rw [getElem!_pos r 0 hr]
        exact Array.getElem_mem_toList hr
      exact ih (divmodRemainder_degreesBound r g
        (Generated.StrictHensel.divmodCoefficient r
          (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int))
        (m : Int) hr hg hrBound hgHead hrDegree hgDegree hactive.2)

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

/-- The remainder returned by the concrete finite division trace is exactly
the state at which the source while-condition became false.  This executable
fact is stronger than the quotient equation alone and is needed to show that
the Hensel correction cannot change a lifted factor's leading term. -/
theorem Generated.StrictHensel.DivmodTrace.output_inactive
    (g : SparsePolyZZ) (inverse m : Int) :
    ∀ {r q : SparsePolyZZ}
      (trace : Generated.StrictHensel.DivmodTrace g inverse m r q),
      let output := Generated.StrictHensel.divmodLoop g inverse m trace
      ¬(!output.2.isEmpty && get_deg output.2 ≥ get_deg g) := by
  intro r q trace
  induction trace with
  | done r q inactive => exact inactive
  | vanished r q active zero next ih => exact ih
  | subtract r q active nonzero next ih => exact ih

/-- Explicit-output form of the actual modular long-division stopping
condition.  No Euclidean-division oracle is used: the result is obtained by
inverting the generated raw call and following its stored finite trace. -/
theorem __upoly_divmod_mod_raw_ir_remainder_inactive_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g q r : SparsePolyZZ) (m : Int)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f g m = .ok (q, r)) :
    ¬(!r.isEmpty && get_deg r ≥ get_deg g) := by
  simp only [Generated.StrictHensel.__upoly_divmod_mod_raw_ir,
    Generated.StrictHensel.__upoly_mod_coeff_raw_ir] at hrun
  split at hrun <;> try contradiction
  next hnotEmpty =>
    split at hrun <;> try contradiction
    next hvalid =>
      split at hrun <;> try contradiction
      next hinvertible =>
        have houtput := Generated.StrictHensel.DivmodTrace.output_inactive
          g (ZZ.invert 0 g[0]!.2 m).2 m
          (termination.trace f g
            (Generated.StrictHensel.modCoeffOutput f m) m rfl hvalid)
        have heq := Except.ok.inj hrun
        rw [heq] at houtput
        exact houtput

/-- The remainder of a successful generated modular-division call preserves
the concrete machine-exponent bound of its input.  This follows the supplied
well-founded trace; it does not replace the generated execution by an
abstract Euclidean remainder. -/
theorem __upoly_divmod_mod_raw_ir_remainder_degreesBound_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g q r : SparsePolyZZ) (m : Nat)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g) (hfBound : DegreesBound f)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f g (m : Int) = .ok (q, r)) :
    DegreesBound r := by
  simp only [Generated.StrictHensel.__upoly_divmod_mod_raw_ir,
    Generated.StrictHensel.__upoly_mod_coeff_raw_ir] at hrun
  split at hrun <;> try contradiction
  next hnotEmpty =>
    split at hrun <;> try contradiction
    next hvalid =>
      split at hrun <;> try contradiction
      next hinvertible =>
        let trace := termination.trace f g
          (Generated.StrictHensel.modCoeffOutput f (m : Int)) (m : Int)
          rfl hvalid
        have hbound := divmodLoop_remainder_degreesBound m g hg hgDegree
          hgHead trace (modCoeffOutput_degreesBound f (m : Int) hfBound)
        have heq := Except.ok.inj hrun
        rw [heq] at hbound
        exact hbound

/-- On the nonempty branch, the concrete stopping condition is the strict
source-degree inequality used by Hensel correction. -/
theorem __upoly_divmod_mod_raw_ir_remainder_get_deg_lt_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g q r : SparsePolyZZ) (m : Int)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f g m = .ok (q, r)) (hr : 0 < r.size) :
    get_deg r < get_deg g := by
  have hinactive :=
    __upoly_divmod_mod_raw_ir_remainder_inactive_of_run termination
      f g q r m hrun
  by_contra hnot
  have hnotInt : ¬(get_deg r).toInt < (get_deg g).toInt := by
    simpa only [Int64.lt_iff_toInt_lt] using hnot
  have hgeInt : (get_deg g).toInt ≤ (get_deg r).toInt :=
    le_of_not_gt hnotInt
  have hge : get_deg r ≥ get_deg g := by
    simpa only [Int64.le_iff_toInt_le] using hgeInt
  apply hinactive
  simp only [Bool.and_eq_true]
  exact ⟨by simp [Array.isEmpty, Nat.ne_of_gt hr], decide_eq_true hge⟩

/-- Every concrete term returned in the generated modular remainder lies
strictly below the divisor's physical head.  This combines the actual loop
stop with canonical sparse ordering and exact signed machine-degree reads. -/
theorem __upoly_divmod_mod_raw_ir_remainder_terms_below_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g q r : SparsePolyZZ) (m : Int)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f g m = .ok (q, r))
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hrCanonical : StrictPolynomialMod.SparsePolyZZCanonical r)
    (hrDegrees : DegreesBound r) :
    ∀ term ∈ r.toList, term.1.deg < g[0]!.1.deg := by
  intro term hterm
  have hr : 0 < r.size := by
    have hlength : 0 < r.toList.length := List.length_pos_iff.mpr
      (List.ne_nil_of_mem hterm)
    simpa using hlength
  have hrHeadMem : r[0]! ∈ r.toList := by
    rw [getElem!_pos r 0 hr]
    exact Array.getElem_mem_toList hr
  have hrDegree : r[0]!.1.deg < 2 ^ 63 := hrDegrees r[0]! hrHeadMem
  have hgetLt :=
    __upoly_divmod_mod_raw_ir_remainder_get_deg_lt_of_run termination
      f g q r m hrun hr
  have hheadLt : r[0]!.1.deg < g[0]!.1.deg := by
    have hgetLtInt : (get_deg r).toInt < (get_deg g).toInt := by
      simpa only [Int64.lt_iff_toInt_lt] using hgetLt
    rw [get_deg_toInt_eq_head r hr hrDegree,
      get_deg_toInt_eq_head g hg hgDegree] at hgetLtInt
    exact_mod_cast hgetLtInt
  exact lt_of_le_of_lt
    (canonical_cursor_degree_le r 0 hr hrCanonical term (by
      simpa using hterm)) hheadLt

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
    (hvalid : termination.inputValid g
      (Generated.StrictHensel.modCoeffOutput f m))
    (hinvert : (ZZ.invert 0 g[0]!.2 m).1 = true) :
    ∃ q r,
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination f g m =
        .ok (q, r) := by
  let reduced := Generated.StrictHensel.modCoeffOutput f m
  let output := Generated.StrictHensel.divmodLoop g
    (ZZ.invert 0 g[0]!.2 m).2 m
      (termination.trace f g reduced m (by rfl) hvalid)
  refine ⟨output.1, output.2, ?_⟩
  simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hg, hvalid, hinvert,
    Generated.StrictHensel.__upoly_mod_coeff_raw_ir, reduced, output]

/-- Full semantic refinement of the generated C++ modular long-division
entry.  The returned concrete arrays satisfy the quotient/remainder equation;
execution is the exact finite trace and not an L2 division call. -/
theorem __upoly_divmod_mod_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g : SparsePolyZZ) (m : Nat)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g) (hfBound : DegreesBound f)
    (hvalid : termination.inputValid g
      (Generated.StrictHensel.modCoeffOutput f (m : Int)))
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true) :
    ∃ q r,
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination f g
          (m : Int) = .ok (q, r) ∧
      toPolyMod m r + toPolyMod m q * toPolyMod m g = toPolyMod m f := by
  let reduced := Generated.StrictHensel.modCoeffOutput f (m : Int)
  let trace := termination.trace f g reduced (m : Int) (by rfl) hvalid
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
      hvalid, hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
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
  have hvalid : termination.inputValid g
      (Generated.StrictHensel.modCoeffOutput f (m : Int)) := by
    simp only [Generated.StrictHensel.__upoly_divmod_mod_raw_ir,
      Generated.StrictHensel.__upoly_mod_coeff_raw_ir] at hrun
    split at hrun <;> try contradiction
    split at hrun <;> try contradiction
    next h => exact h
  rcases __upoly_divmod_mod_raw_ir_refines termination f g m hg hgDegree
      hgHead hfBound hvalid hinvert with
    ⟨actualQ, actualR, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

/-- A successful generated modular-divmod call returns the remainder of its
actual finite trace.  Canonicality follows from that trace, not from a
specification-level quotient/remainder oracle. -/
theorem __upoly_divmod_mod_raw_ir_remainder_canonical_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g q r : SparsePolyZZ) (m : Int)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f g m = .ok (q, r)) :
    StrictPolynomialMod.SparsePolyZZCanonical r := by
  have hgFalse : g.isEmpty = false := by
    simp only [Generated.StrictHensel.__upoly_divmod_mod_raw_ir,
      Generated.StrictHensel.__upoly_mod_coeff_raw_ir] at hrun
    split at hrun
    · contradiction
    · rename_i hnotEmpty
      exact Bool.eq_false_of_not_eq_true hnotEmpty
  have hvalid : termination.inputValid g
      (Generated.StrictHensel.modCoeffOutput f m) := by
    simp only [Generated.StrictHensel.__upoly_divmod_mod_raw_ir,
      Generated.StrictHensel.__upoly_mod_coeff_raw_ir] at hrun
    split at hrun <;> try contradiction
    split at hrun <;> try contradiction
    next h => exact h
  have hinvert : (ZZ.invert 0 g[0]!.2 m).1 = true := by
    by_contra hinvert
    have hinvertFalse : (ZZ.invert 0 g[0]!.2 m).1 = false :=
      Bool.eq_false_of_not_eq_true hinvert
    simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hgFalse,
      hvalid, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      hinvertFalse] at hrun
  let reduced := Generated.StrictHensel.modCoeffOutput f m
  let trace := termination.trace f g reduced m (by rfl) hvalid
  let actual := Generated.StrictHensel.divmodLoop g
    (ZZ.invert 0 g[0]!.2 m).2 m trace
  have hactual :
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination f g m =
        .ok (actual.1, actual.2) := by
    simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hgFalse,
      hvalid, hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      reduced, trace, actual]
  have hcanonical : StrictPolynomialMod.SparsePolyZZCanonical actual.2 := by
    exact Refinement.StrictHensel.Generated.StrictHensel.DivmodTrace.remainder_canonical
      trace
      (modCoeffOutput_canonical f m hf) hg
  rw [hrun] at hactual
  cases hactual
  exact hcanonical

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

/-- The two factors produced by the first concrete Hensel phase remain strict
sparse representations.  Every intermediate below is fixed by a successful
generated raw operation; in particular the modular remainder is obtained from
the actual well-founded division trace. -/
theorem henselFactorCorrection_canonical_from_raw_runs
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat) (hm : 0 < m)
    (gh difference se q r te qg tauRaw tau gRaw gNew hRaw hNew :
      SparsePolyZZ)
    (hfCanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hgCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.g)
    (hhCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.h)
    (hgh : (strictHenselRawOps termination).mul node.g node.h = .ok gh)
    (hdifference :
      (strictHenselRawOps termination).sub f gh = .ok difference)
    (hse : (strictHenselRawOps termination).mul node.s
      (Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)) =
        .ok se)
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
      (m ^ 2 : Int) = .ok hNew) :
    StrictPolynomialMod.SparsePolyZZCanonical gNew ∧
      StrictPolynomialMod.SparsePolyZZCanonical hNew := by
  have hmInt : (m : Int) ≠ 0 := by omega
  have hghCanonical := strictHenselRawOps_mul_canonical_of_run
    termination node.g node.h gh hgh
  have hdifferenceCanonical := strictHenselRawOps_sub_canonical_of_run
    termination f gh difference hfCanonical hghCanonical hdifference
  have heCanonical := divideThenReduceCoeffs_canonical difference (m : Int)
    hdifferenceCanonical
  have hseCanonical := strictHenselRawOps_mul_canonical_of_run
    termination node.s _ se hse
  have hrCanonical := __upoly_divmod_mod_raw_ir_remainder_canonical_of_run
    termination se node.h q r (m : Int) hseCanonical hhCanonical hdivmod
  have hteCanonical := strictHenselRawOps_mul_canonical_of_run
    termination node.t _ te hte
  have hqgCanonical := strictHenselRawOps_mul_canonical_of_run
    termination q node.g qg hqg
  have htauRawCanonical := strictHenselRawOps_add_canonical_of_run
    termination te qg tauRaw hteCanonical hqgCanonical htauRaw
  have htauCanonical := modCoeffOutput_canonical_of_run tauRaw tau (m : Int)
    htauRawCanonical htau
  have hscaledTauCanonical := scaleCoeffs_canonical tau (m : Int)
    htauCanonical hmInt
  have hgRawCanonical := strictHenselRawOps_add_canonical_of_run
    termination node.g _ gRaw hgCanonical hscaledTauCanonical hgRaw
  have hgNewCanonical := modCoeffOutput_canonical_of_run gRaw gNew
    (m ^ 2 : Int) hgRawCanonical hgNew
  have hscaledRCanonical := scaleCoeffs_canonical r (m : Int)
    hrCanonical hmInt
  have hhRawCanonical := strictHenselRawOps_add_canonical_of_run
    termination node.h _ hRaw hhCanonical hscaledRCanonical hhRaw
  have hhNewCanonical := modCoeffOutput_canonical_of_run hRaw hNew
    (m ^ 2 : Int) hhRawCanonical hhNew
  exact ⟨hgNewCanonical, hhNewCanonical⟩

/-- The physical coefficient-one head of the right factor survives one
successful generated Hensel factor phase.  The proof follows the actual
modular division remainder and the actual scale/merge/reduce calls. -/
theorem henselFactorCorrection_preserves_h_one_head_from_raw_runs
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (m : Nat) (hm : 1 < m)
    (e se q r hRaw hNew : SparsePolyZZ) (head : UMonomial)
    (tail : List (UMonomial × Int))
    (hhList : node.h.toList = (head, 1) :: tail)
    (hhCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.h)
    (hh : 0 < node.h.size) (hhDegree : node.h[0]!.1.deg < 2 ^ 63)
    (hhHead : HeadDominates node.h) (hseBound : DegreesBound se)
    (hseRun : (strictHenselRawOps termination).mul node.s
      e = .ok se)
    (hdivmod : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      se node.h (m : Int) = .ok (q, r))
    (hhRaw : (strictHenselRawOps termination).add node.h
      (Generated.StrictHensel.scaleCoeffs r (m : Int)) = .ok hRaw)
    (hhNew : Generated.StrictHensel.__upoly_mod_coeff_raw_ir hRaw
      (m ^ 2 : Int) = .ok hNew) :
    ∃ suffix, hNew.toList = (head, 1) :: suffix := by
  have hseCanonical := strictHenselRawOps_mul_canonical_of_run
    termination node.s e se hseRun
  have hrCanonical := __upoly_divmod_mod_raw_ir_remainder_canonical_of_run
    termination se node.h q r (m : Int) hseCanonical hhCanonical hdivmod
  have hrDegrees :=
    __upoly_divmod_mod_raw_ir_remainder_degreesBound_of_run termination
      se node.h q r m hh hhDegree hhHead hseBound hdivmod
  have hrBelow := __upoly_divmod_mod_raw_ir_remainder_terms_below_of_run
    termination se node.h q r (m : Int) hdivmod hh hhDegree hrCanonical
      hrDegrees
  have hhHeadEq : node.h[0]!.1.deg = head.deg := by
    have hhead := congrArg List.head? hhList
    have hget : node.h[0]? = some (head, 1) := by
      rw [List.head?_eq_getElem?, Array.getElem?_toList] at hhead
      exact hhead
    have hzero : node.h[0]! = (head, 1) := by
      rw [getElem!_pos node.h 0 hh]
      exact Option.some.inj ((Array.getElem?_eq_getElem hh).symm.trans hget)
    simpa [hzero]
  have hrBelowHead : ∀ term ∈ r.toList, term.1.deg < head.deg := by
    simpa [hhHeadEq] using hrBelow
  have hphysical := henselHCorrection_preserves_one_head node.h r m head tail
    hhList hrBelowHead hm
  have hhRawEq : hRaw = Generated.StrictHensel.pairVecAddLoop node.h
      (Generated.StrictHensel.scaleCoeffs r (m : Int)) 0 0 #[] := by
    simpa only [strictHenselRawOps,
      Generated.StrictHensel.__upoly_add_raw_ir] using
        (Except.ok.inj hhRaw).symm
  have hhNewEq : hNew = Generated.StrictHensel.modCoeffOutput hRaw
      (m ^ 2 : Int) := by
    simpa only [Generated.StrictHensel.__upoly_mod_coeff_raw_ir] using
      (Except.ok.inj hhNew).symm
  rw [hhNewEq, hhRawEq]
  exact hphysical

set_option maxHeartbeats 0 in
/-- The complete generated factor phase refines its L2 factor invariant.
All local values below are definitionally the outputs of the concrete
generated operations, including the result of the supplied well-founded
division trace. -/
theorem __hensel_step_factor_phase_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat) (hm : 0 < m)
    (hfCanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hgCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.g)
    (hhCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.h)
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
    (hdivmodValid :
      let gh := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
      let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
      let e := Generated.StrictHensel.divideThenReduceCoeffs
        difference (m : Int)
      let se := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
      termination.inputValid node.h
        (Generated.StrictHensel.modCoeffOutput se (m : Int)))
    (hinvariant : HenselNodeInvariant f m node) :
    ∃ factorNode,
      Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
          (strictHenselRawOps termination) node f (m : Int) = .ok factorNode ∧
      toPolyMod (m ^ 2) factorNode.g *
          toPolyMod (m ^ 2) factorNode.h = toPolyMod (m ^ 2) f ∧
      toPolyMod m factorNode.g = toPolyMod m node.g ∧
      toPolyMod m factorNode.h = toPolyMod m node.h ∧
      factorNode.s = node.s ∧ factorNode.t = node.t ∧
      StrictPolynomialMod.SparsePolyZZCanonical factorNode.g ∧
      StrictPolynomialMod.SparsePolyZZCanonical factorNode.h ∧
      ∀ head tail, node.h.toList = (head, 1) :: tail → 1 < m →
        ∃ suffix, factorNode.h.toList = (head, 1) :: suffix := by
  let gh := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
  let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
  let e := Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)
  let se := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
  let reduced := Generated.StrictHensel.modCoeffOutput se (m : Int)
  let trace := termination.trace se node.h reduced (m : Int) (by rfl)
    hdivmodValid
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
      hdivmodValid, hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      reduced, trace, qr]
    rw [dif_pos hdivmodValid]
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
  have hcanonical := henselFactorCorrection_canonical_from_raw_runs
    termination node f m hm gh difference se qr.1 qr.2 te qg tauRaw tau gRaw
    gNew hRaw hNew hfCanonical hgCanonical hhCanonical hghRun hdifferenceRun
    hseRun hdivmodRun hteRun hqgRun htauRawRun htauRun hgRawRun hgNewRun
    hhRawRun hhNewRun
  have hhPhysical : ∀ head tail,
      node.h.toList = (head, 1) :: tail → 1 < m →
      ∃ suffix, factorNode.h.toList = (head, 1) :: suffix := by
    intro head tail hhList hmOne
    exact henselFactorCorrection_preserves_h_one_head_from_raw_runs
      termination node m hmOne e se qr.1 qr.2 hRaw hNew head tail hhList
      hhCanonical hh hhDegree hhHead
      (by simpa [gh, difference, e, se] using hseBound) hseRun hdivmodRun
      hhRawRun hhNewRun
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
  · exact ⟨hsemantic.1, hsemantic.2.1, hsemantic.2.2, rfl, rfl,
      hcanonical.1, hcanonical.2, hhPhysical⟩

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
    (hdivmodValid :
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
      termination.inputValid factorNode.h
        (Generated.StrictHensel.modCoeffOutput sep (m : Int)))
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
    hdivmodValid
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
      hdivmodValid, hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      reduced, trace, qr]
    rw [dif_pos hdivmodValid]
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
  divmodValid :
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
    termination.inputValid factorNode.h
      (Generated.StrictHensel.modCoeffOutput sep (m : Int))

/-- Complete algorithm invariant for the generated quadratic Hensel step.
It records only source safety, representation bounds, exact divisibility, and
the input L2 invariant.  The intermediate factor node is universally tied to
the unique generated first-phase execution. -/
structure HenselStepRefinementInvariant
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat) : Prop where
  positiveModulus : 0 < m
  inputInvariant : HenselNodeInvariant f m node
  targetCanonical : StrictPolynomialMod.SparsePolyZZCanonical f
  gCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.g
  hCanonical : StrictPolynomialMod.SparsePolyZZCanonical node.h
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
  factorDivmodValid :
    let gh := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
    let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
    let e := Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)
    let se := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
    termination.inputValid node.h
      (Generated.StrictHensel.modCoeffOutput se (m : Int))
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
      hinvariant.positiveModulus hinvariant.targetCanonical
      hinvariant.gCanonical hinvariant.hCanonical
      hinvariant.hNonempty hinvariant.hDegree
      hinvariant.hHead hinvariant.hInvertible
      hinvariant.factorDifferenceDivisible hinvariant.factorSeBound
      hinvariant.factorDivmodValid
      hinvariant.inputInvariant with
    ⟨factorNode, hfactorRun, hfactorProduct, hgPreserved, hhPreserved,
      hsUnchanged, htUnchanged, hgCanonical, hhCanonical, hhPhysical⟩
  have hready := hinvariant.bezoutReady factorNode hfactorRun
  have hfactorBezout :
      toPolyMod m factorNode.s * toPolyMod m factorNode.g +
        toPolyMod m factorNode.t * toPolyMod m factorNode.h = 1 := by
    rw [hsUnchanged, htUnchanged, hgPreserved, hhPreserved]
    exact hinvariant.inputInvariant.2
  rcases __hensel_step_bezout_phase_raw_ir_refines termination factorNode m
      hinvariant.positiveModulus hready.hNonempty hready.hDegree hready.hHead
      hready.hInvertible hready.differenceDivisible hready.sepBound
      hready.divmodValid
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
    constructor
    · calc
        toPolyMod m output.t = toPolyMod m factorNode.t := htPreserved
        _ = toPolyMod m node.t := by rw [htUnchanged]
    constructor
    · rw [hgUnchanged]
      exact hgCanonical
    constructor
    · rw [hhUnchanged]
      exact hhCanonical
    · intro head tail hhList hmOne
      rw [hhUnchanged]
      exact hhPhysical head tail hhList hmOne

/-- The first contiguous source phase has a genuine raw-to-safe execution
bridge.  Its only possible source assertion is the modular division by `h`;
all generated heap arithmetic and coefficient-reduction calls are total. -/
theorem __hensel_step_factor_phase_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : ZZ)
    (hh : node.h.isEmpty = false)
    (hvalid :
      let gh := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
      let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
      let e := Generated.StrictHensel.divideThenReduceCoeffs difference m
      let se := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
      termination.inputValid node.h
        (Generated.StrictHensel.modCoeffOutput se m))
    (hinvert : (ZZ.invert 0 node.h[0]!.2 m).1 = true) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
          (strictHenselRawOps termination) node f m = .ok output := by
  let gh := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
  let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
  let e := Generated.StrictHensel.divideThenReduceCoeffs difference m
  let se := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
  rcases __upoly_divmod_mod_raw_ir_terminates termination se node.h m hh
      hvalid hinvert with ⟨q, r, hdivmod⟩
  unfold Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
    strictHenselRawOps Generated.StrictHensel.__upoly_mul_raw_ir
    Generated.StrictHensel.__upoly_add_raw_ir
    Generated.StrictHensel.__upoly_sub_raw_ir
  simp only [Bind.bind, Except.bind]
  dsimp [se, e, difference, gh] at hdivmod
  rw [hdivmod]
  simp only [Bind.bind, Except.bind]
  exact ⟨_, rfl⟩

/-- Matching raw-to-safe bridge for the second contiguous source phase. -/
theorem __hensel_step_bezout_phase_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (factorNode : HenselNode) (m : ZZ)
    (hh : factorNode.h.isEmpty = false)
    (hvalid :
      let sg := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g) #[]
      let th := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h) #[]
      let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
        (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
      let difference := Generated.StrictHensel.pairVecSubLoop oneMinusSg th 0 0 #[]
      let ep := Generated.StrictHensel.divideThenReduceCoeffs difference m
      let sep := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.s ep) #[]
      termination.inputValid factorNode.h
        (Generated.StrictHensel.modCoeffOutput sep m))
    (hinvert : (ZZ.invert 0 factorNode.h[0]!.2 m).1 = true) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir
          (strictHenselRawOps termination) factorNode m = .ok output := by
  let sg := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g) #[]
  let th := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h) #[]
  let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
    (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
  let difference := Generated.StrictHensel.pairVecSubLoop oneMinusSg th 0 0 #[]
  let ep := Generated.StrictHensel.divideThenReduceCoeffs difference m
  let sep := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.s ep) #[]
  rcases __upoly_divmod_mod_raw_ir_terminates termination sep factorNode.h m
      hh hvalid hinvert with ⟨q, r, hdivmod⟩
  unfold Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir
    strictHenselRawOps Generated.StrictHensel.__upoly_mul_raw_ir
    Generated.StrictHensel.__upoly_add_raw_ir
    Generated.StrictHensel.__upoly_sub_raw_ir
  simp only [Bind.bind, Except.bind]
  dsimp [sep, ep, difference, oneMinusSg, th, sg] at hdivmod
  rw [hdivmod]
  simp only [Bind.bind, Except.bind]
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
  factorDivmodValid :
    let gh := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
    let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
    let e := Generated.StrictHensel.divideThenReduceCoeffs difference m
    let se := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
    termination.inputValid node.h
      (Generated.StrictHensel.modCoeffOutput se m)
  factorHReady : ∀ factorNode,
    Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
        (strictHenselRawOps termination) node f m = .ok factorNode →
    factorNode.h.isEmpty = false ∧
      (ZZ.invert 0 factorNode.h[0]!.2 m).1 = true ∧
      (let sg := Generated.StrictHensel.pairVecMulHeapLoop
          (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g) #[]
       let th := Generated.StrictHensel.pairVecMulHeapLoop
          (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h) #[]
       let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
          (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
       let difference := Generated.StrictHensel.pairVecSubLoop oneMinusSg th 0 0 #[]
       let ep := Generated.StrictHensel.divideThenReduceCoeffs difference m
       let sep := Generated.StrictHensel.pairVecMulHeapLoop
          (Generated.StrictHensel.pairVecMulProducts factorNode.s ep) #[]
       termination.inputValid factorNode.h
          (Generated.StrictHensel.modCoeffOutput sep m))

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
      hinvariant.inputHNonempty hinvariant.factorDivmodValid
      hinvariant.inputHInvertible with
    ⟨factorNode, hfactor⟩
  have hready := hinvariant.factorHReady factorNode hfactor
  rcases __hensel_step_bezout_phase_raw_ir_terminates termination factorNode m
      hready.1 hready.2.2 hready.2.1 with ⟨output, hbezout⟩
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

/-- The exact generated Hensel step updates only `g`, `h`, `s`, and `t`.
Consequently every successful raw execution preserves the two source child
indices and leaf interval fields. -/
theorem __hensel_step_raw_ir_preserves_topology
    (ops : Generated.StrictHensel.HenselStepRawOps)
    (node output : HenselNode) (f : SparsePolyZZ) (m : ZZ)
    (hrun : Generated.StrictHensel.__hensel_step_raw_ir ops node f m =
      .ok output) :
    output.left = node.left ∧ output.right = node.right ∧
      output.leaf_start = node.leaf_start ∧ output.leaf_end = node.leaf_end := by
  have bind_ok {α β : Type} {value : Except RawFault α}
      {next : α → Except RawFault β} {result : β}
      (h : value >>= next = .ok result) :
      ∃ intermediate, value = .ok intermediate ∧
        next intermediate = .ok result := by
    cases value with
    | error fault =>
        change (.error fault : Except RawFault β) = .ok result at h
        contradiction
    | ok intermediate => exact ⟨intermediate, rfl, h⟩
  rw [Generated.StrictHensel.__hensel_step_raw_ir] at hrun
  rcases bind_ok hrun with ⟨factorNode, hfactor, hbezout⟩
  have hfactorTopology :
      factorNode.left = node.left ∧ factorNode.right = node.right ∧
      factorNode.leaf_start = node.leaf_start ∧
      factorNode.leaf_end = node.leaf_end := by
    simp only [Generated.StrictHensel.__hensel_step_factor_phase_raw_ir,
      bind, Except.bind] at hfactor
    repeat' split at hfactor <;> simp_all
    injection hfactor with hfactorEq
    subst factorNode
    simp
  have hbezoutTopology :
      output.left = factorNode.left ∧ output.right = factorNode.right ∧
      output.leaf_start = factorNode.leaf_start ∧
      output.leaf_end = factorNode.leaf_end := by
    simp only [Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir,
      bind, Except.bind] at hbezout
    repeat' split at hbezout <;> simp_all
    injection hbezout with hbezoutEq
    subst output
    simp
  exact ⟨hbezoutTopology.1.trans hfactorTopology.1,
    hbezoutTopology.2.1.trans hfactorTopology.2.1,
    hbezoutTopology.2.2.1.trans hfactorTopology.2.2.1,
    hbezoutTopology.2.2.2.trans hfactorTopology.2.2.2⟩

def HenselNodeTopologyEq (left right : HenselNode) : Prop :=
  left.left = right.left ∧ left.right = right.right ∧
    left.leaf_start = right.leaf_start ∧ left.leaf_end = right.leaf_end

def HenselArrayTopologyEq
    (left right : Array HenselNode) : Prop :=
  left.size = right.size ∧
    ∀ index (hleft : index < left.size) (hright : index < right.size),
      HenselNodeTopologyEq left[index] right[index]

theorem HenselNodeTopologyEq.refl (node : HenselNode) :
    HenselNodeTopologyEq node node := by
  simp [HenselNodeTopologyEq]

theorem HenselNodeTopologyEq.of_eq {left right : HenselNode}
    (h : left = right) : HenselNodeTopologyEq left right := by
  subst right
  exact .refl _

theorem HenselNodeTopologyEq.trans {first second third : HenselNode}
    (hfirst : HenselNodeTopologyEq first second)
    (hsecond : HenselNodeTopologyEq second third) :
    HenselNodeTopologyEq first third := by
  simp only [HenselNodeTopologyEq] at hfirst hsecond ⊢
  exact ⟨hfirst.1.trans hsecond.1,
    hfirst.2.1.trans hsecond.2.1,
    hfirst.2.2.1.trans hsecond.2.2.1,
    hfirst.2.2.2.trans hsecond.2.2.2⟩

theorem HenselNodeTopologyEq.symm {left right : HenselNode}
    (h : HenselNodeTopologyEq left right) :
    HenselNodeTopologyEq right left := by
  simp only [HenselNodeTopologyEq] at h ⊢
  exact ⟨h.1.symm, h.2.1.symm, h.2.2.1.symm, h.2.2.2.symm⟩

theorem HenselArrayTopologyEq.refl (nodes : Array HenselNode) :
    HenselArrayTopologyEq nodes nodes := by
  exact ⟨rfl, fun index hleft hright => .refl _⟩

theorem HenselArrayTopologyEq.trans {first second third : Array HenselNode}
    (hfirst : HenselArrayTopologyEq first second)
    (hsecond : HenselArrayTopologyEq second third) :
    HenselArrayTopologyEq first third := by
  refine ⟨hfirst.1.trans hsecond.1, ?_⟩
  intro index hfirstIndex hthirdIndex
  have hsecondIndex : index < second.size := by simpa [← hfirst.1]
  exact (hfirst.2 index hfirstIndex hsecondIndex).trans
    (hsecond.2 index hsecondIndex hthirdIndex)

theorem henselArrayTopologyEq_set
    (nodes : Array HenselNode) (index : Nat) (value : HenselNode)
    (hindex : index < nodes.size)
    (hvalue : HenselNodeTopologyEq nodes[index]! value) :
    HenselArrayTopologyEq nodes (nodes.set! index value) := by
  refine ⟨by simp, ?_⟩
  intro other hbefore hafter
  by_cases heq : other = index
  · subst other
    simpa [Array.getElem_set, hindex] using hvalue
  · simpa only [Array.set!_eq_setIfInBounds, Array.setIfInBounds_def,
      dif_pos hindex, Array.getElem_set_ne hindex hbefore (Ne.symm heq)] using
      HenselNodeTopologyEq.refl nodes[other]

/-- Algebraic fields of `after` reduce to the corresponding fields of
`before` modulo `m`.  Pointer/topology fields are deliberately excluded: they
are tracked independently by `HenselNodeTopologyEq`. -/
def HenselNodeReduces (m : Nat) (before after : HenselNode) : Prop :=
  toPolyMod m after.g = toPolyMod m before.g ∧
  toPolyMod m after.h = toPolyMod m before.h ∧
  toPolyMod m after.s = toPolyMod m before.s ∧
  toPolyMod m after.t = toPolyMod m before.t

def HenselArrayReduces (m : Nat)
    (before after : Array HenselNode) : Prop :=
  before.size = after.size ∧
  ∀ index (hbefore : index < before.size) (hafter : index < after.size),
    HenselNodeReduces m before[index] after[index]

/-- Representation invariant for the two factor fields stored at a concrete
C++ Hensel node. -/
def HenselNodeCanonical (node : HenselNode) : Prop :=
  StrictPolynomialMod.SparsePolyZZCanonical node.g ∧
  StrictPolynomialMod.SparsePolyZZCanonical node.h

/-- Every allocated node in the concrete Hensel array stores canonical sparse
factor fields. -/
def HenselArrayCanonical (nodes : Array HenselNode) : Prop :=
  ∀ index (hindex : index < nodes.size), HenselNodeCanonical nodes[index]

/-- A physical sparse factor is monic when coefficient one is stored at its
first array position.  This representation-level property is intentionally
stronger than mathematical monicity of its decoded polynomial. -/
def HasPhysicalOneHead (factor : SparsePolyZZ) : Prop :=
  ∃ head tail, factor.toList = (head, 1) :: tail

def HenselArrayHOneHead (nodes : Array HenselNode) : Prop :=
  ∀ index (hindex : index < nodes.size), HasPhysicalOneHead nodes[index].h

theorem henselArrayHOneHead_set
    (nodes : Array HenselNode) (index : Nat) (value : HenselNode)
    (hindex : index < nodes.size) (hnodes : HenselArrayHOneHead nodes)
    (hvalue : HasPhysicalOneHead value.h) :
    HenselArrayHOneHead (nodes.set! index value) := by
  intro other hother
  have hotherBefore : other < nodes.size := by simpa using hother
  by_cases heq : other = index
  · subst other
    simpa [Array.getElem_set, hindex] using hvalue
  · simpa only [Array.set!_eq_setIfInBounds, Array.setIfInBounds_def,
      dif_pos hindex,
      Array.getElem_set_ne hindex hotherBefore (Ne.symm heq)] using
      hnodes other hotherBefore

theorem henselArrayCanonical_set
    (nodes : Array HenselNode) (index : Nat) (value : HenselNode)
    (hindex : index < nodes.size) (hnodes : HenselArrayCanonical nodes)
    (hvalue : HenselNodeCanonical value) :
    HenselArrayCanonical (nodes.set! index value) := by
  intro other hother
  have hotherBefore : other < nodes.size := by simpa using hother
  by_cases heq : other = index
  · subst other
    simpa [Array.getElem_set, hindex] using hvalue
  · simpa only [Array.set!_eq_setIfInBounds, Array.setIfInBounds_def,
      dif_pos hindex,
      Array.getElem_set_ne hindex hotherBefore (Ne.symm heq)] using
      hnodes other hotherBefore

theorem HenselNodeReduces.refl (m : Nat) (node : HenselNode) :
    HenselNodeReduces m node node := by
  simp [HenselNodeReduces]

theorem HenselNodeReduces.trans {m : Nat} {first second third : HenselNode}
    (hfirst : HenselNodeReduces m first second)
    (hsecond : HenselNodeReduces m second third) :
    HenselNodeReduces m first third := by
  simp only [HenselNodeReduces] at hfirst hsecond ⊢
  exact ⟨hsecond.1.trans hfirst.1,
    hsecond.2.1.trans hfirst.2.1,
    hsecond.2.2.1.trans hfirst.2.2.1,
    hsecond.2.2.2.trans hfirst.2.2.2⟩

theorem HenselArrayReduces.refl (m : Nat) (nodes : Array HenselNode) :
    HenselArrayReduces m nodes nodes := by
  exact ⟨rfl, fun _ _ _ => HenselNodeReduces.refl m _⟩

theorem HenselArrayReduces.trans {m : Nat}
    {first second third : Array HenselNode}
    (hfirst : HenselArrayReduces m first second)
    (hsecond : HenselArrayReduces m second third) :
    HenselArrayReduces m first third := by
  refine ⟨hfirst.1.trans hsecond.1, ?_⟩
  intro index hfirstIndex hthirdIndex
  have hsecondIndex : index < second.size := by simpa [← hfirst.1]
  exact (hfirst.2 index hfirstIndex hsecondIndex).trans
    (hsecond.2 index hsecondIndex hthirdIndex)

theorem HenselNodeReduces.of_dvd {smaller larger : Nat}
    (hdiv : smaller ∣ larger) {before after : HenselNode}
    (hreduce : HenselNodeReduces larger before after) :
    HenselNodeReduces smaller before after := by
  simp only [HenselNodeReduces] at hreduce ⊢
  have project (left right : SparsePolyZZ)
      (h : toPolyMod larger left = toPolyMod larger right) :
      toPolyMod smaller left = toPolyMod smaller right := by
    have hmapped := congrArg
      (Polynomial.map (ZMod.castHom hdiv (ZMod smaller))) h
    simpa only [map_toPolyMod_of_dvd hdiv] using hmapped
  exact ⟨project after.g before.g hreduce.1,
    project after.h before.h hreduce.2.1,
    project after.s before.s hreduce.2.2.1,
    project after.t before.t hreduce.2.2.2⟩

theorem HenselArrayReduces.of_dvd {smaller larger : Nat}
    (hdiv : smaller ∣ larger) {before after : Array HenselNode}
    (hreduce : HenselArrayReduces larger before after) :
    HenselArrayReduces smaller before after := by
  refine ⟨hreduce.1, ?_⟩
  intro index hbefore hafter
  exact (hreduce.2 index hbefore hafter).of_dvd hdiv

theorem henselArrayReduces_set
    (m : Nat) (nodes : Array HenselNode) (index : Nat) (value : HenselNode)
    (hindex : index < nodes.size)
    (hvalue : HenselNodeReduces m nodes[index]! value) :
    HenselArrayReduces m nodes (nodes.set! index value) := by
  refine ⟨by simp, ?_⟩
  intro other hbefore hafter
  by_cases heq : other = index
  · subst other
    simpa [Array.getElem_set, hindex] using hvalue
  · simpa only [Array.set!_eq_setIfInBounds, Array.setIfInBounds_def,
      dif_pos hindex, Array.getElem_set_ne hindex hbefore (Ne.symm heq)] using
      HenselNodeReduces.refl m nodes[other]

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

/-- Concrete array indices visited by a finite generated Hensel topology, in
preorder.  This contains no polynomial semantics; it is used to state exact
mutable-array frame properties of the real recursive execution. -/
def henselLiftTreeIndices : Generated.StrictHensel.HenselLiftTree → List Nat
  | .node index none none => [index]
  | .node index (some left) none => index :: henselLiftTreeIndices left
  | .node index none (some right) => index :: henselLiftTreeIndices right
  | .node index (some left) (some right) =>
      index :: henselLiftTreeIndices left ++ henselLiftTreeIndices right

private theorem array_getElem?_setBang_ne {α : Type*}
    (items : Array α) (source target : Nat) (value : α)
    (hsource : source < items.size) (hne : target ≠ source) :
    (items.set! source value)[target]? = items[target]? := by
  exact Array.getElem?_setIfInBounds_ne hne.symm

/-- A successful generated recursive Hensel traversal leaves every array
cell whose index is absent from its concrete finite topology unchanged. -/
theorem HenselLiftRecursiveCorrect.preserves_not_mem
    {termination : Generated.StrictHensel.DivmodTermination} {m : Nat}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes output : Array HenselNode} {target : SparsePolyZZ}
    (hcorrect : HenselLiftRecursiveCorrect termination m tree nodes target
      output) :
    ∀ index, index ∉ henselLiftTreeIndices tree →
      output[index]? = nodes[index]? := by
  induction hcorrect with
  | leaf root nodes stored target inputNode lifted parent hnode hstep
      hstepCorrect hstored hparent =>
      subst stored
      intro index houtside
      have hroot : root < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hne : index ≠ root := by
        simpa [henselLiftTreeIndices] using houtside
      exact array_getElem?_setBang_ne nodes root index lifted hroot hne
  | left root left nodes stored nodesAfterLeft target inputNode lifted parent
      hnode hstep hstepCorrect hstored hleftRun hleftCorrect hparent ih =>
      subst stored
      intro index houtside
      have hroot : root < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hne : index ≠ root := by
        intro heq
        subst index
        exact houtside (by simp [henselLiftTreeIndices])
      have houtsideLeft : index ∉ henselLiftTreeIndices left := by
        have hparts : index ≠ root ∧
            index ∉ henselLiftTreeIndices left := by
          simpa [henselLiftTreeIndices] using houtside
        exact hparts.2
      calc
        nodesAfterLeft[index]? = (nodes.set! root lifted)[index]? :=
          ih index houtsideLeft
        _ = nodes[index]? :=
          array_getElem?_setBang_ne nodes root index lifted hroot hne
  | right root right nodes stored output target inputNode lifted parent hnode
      hstep hstepCorrect hstored hparent hrightRun hrightCorrect ih =>
      subst stored
      intro index houtside
      have hroot : root < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hne : index ≠ root := by
        intro heq
        subst index
        exact houtside (by simp [henselLiftTreeIndices])
      have houtsideRight : index ∉ henselLiftTreeIndices right := by
        have hparts : index ≠ root ∧
            index ∉ henselLiftTreeIndices right := by
          simpa [henselLiftTreeIndices] using houtside
        exact hparts.2
      calc
        output[index]? = (nodes.set! root lifted)[index]? :=
          ih index houtsideRight
        _ = nodes[index]? :=
          array_getElem?_setBang_ne nodes root index lifted hroot hne
  | branch root left right nodes stored nodesAfterLeft output target inputNode
      lifted parent hnode hstep hstepCorrect hstored hleftRun hleftCorrect
      hparent hrightRun hrightCorrect leftIH rightIH =>
      subst stored
      intro index houtside
      have hroot : root < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hne : index ≠ root := by
        intro heq
        subst index
        exact houtside (by simp [henselLiftTreeIndices])
      have houtsideLeft : index ∉ henselLiftTreeIndices left := by
        have hparts : index ≠ root ∧
            index ∉ henselLiftTreeIndices left ∧
              index ∉ henselLiftTreeIndices right := by
          simpa [henselLiftTreeIndices] using houtside
        exact hparts.2.1
      have houtsideRight : index ∉ henselLiftTreeIndices right := by
        have hparts : index ≠ root ∧
            index ∉ henselLiftTreeIndices left ∧
              index ∉ henselLiftTreeIndices right := by
          simpa [henselLiftTreeIndices] using houtside
        exact hparts.2.2
      calc
        output[index]? = nodesAfterLeft[index]? :=
          rightIH index houtsideRight
        _ = (nodes.set! root lifted)[index]? :=
          leftIH index houtsideLeft
        _ = nodes[index]? :=
          array_getElem?_setBang_ne nodes root index lifted hroot hne

theorem HenselLiftRecursiveCorrect.topologyEq
    {termination : Generated.StrictHensel.DivmodTermination} {m : Nat}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes output : Array HenselNode} {target : SparsePolyZZ}
    (hcorrect : HenselLiftRecursiveCorrect termination m tree nodes target
      output) :
    HenselArrayTopologyEq nodes output := by
  induction hcorrect with
  | leaf index nodes stored target inputNode lifted parent hnode hstep
      hstepCorrect hstored hparent =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode := by
          exact Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans
            hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      have htopology := __hensel_step_raw_ir_preserves_topology
        (strictHenselRawOps termination) inputNode lifted target (m : Int)
        hstep
      apply henselArrayTopologyEq_set nodes index lifted hindex
      rw [hinput]
      exact (show HenselNodeTopologyEq lifted inputNode from htopology).symm
  | left index left nodes stored nodesAfterLeft target inputNode lifted parent
      hnode hstep hstepCorrect hstored hleftRun hleftCorrect hparent ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode := by
          exact Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans
            hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      have htopology := __hensel_step_raw_ir_preserves_topology
        (strictHenselRawOps termination) inputNode lifted target (m : Int)
        hstep
      exact (henselArrayTopologyEq_set nodes index lifted hindex (by
        rw [hinput]
        exact (show HenselNodeTopologyEq lifted inputNode from htopology).symm)).trans ih
  | right index right nodes stored output target inputNode lifted parent hnode
      hstep hstepCorrect hstored hparent hrightRun hrightCorrect ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode := by
          exact Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans
            hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      have htopology := __hensel_step_raw_ir_preserves_topology
        (strictHenselRawOps termination) inputNode lifted target (m : Int)
        hstep
      exact (henselArrayTopologyEq_set nodes index lifted hindex (by
        rw [hinput]
        exact (show HenselNodeTopologyEq lifted inputNode from htopology).symm)).trans ih
  | branch index left right nodes stored nodesAfterLeft output target inputNode
      lifted parent hnode hstep hstepCorrect hstored hleftRun hleftCorrect
      hparent hrightRun hrightCorrect leftIH rightIH =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode := by
          exact Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans
            hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      have htopology := __hensel_step_raw_ir_preserves_topology
        (strictHenselRawOps termination) inputNode lifted target (m : Int)
        hstep
      exact ((henselArrayTopologyEq_set nodes index lifted hindex (by
        rw [hinput]
        exact (show HenselNodeTopologyEq lifted inputNode from htopology).symm)).trans
          leftIH).trans rightIH

/-- Every algebraic node field produced by the concrete recursive tree walk
reduces modulo the input modulus to the field stored before that walk. -/
theorem HenselLiftRecursiveCorrect.arrayReduces
    {termination : Generated.StrictHensel.DivmodTermination} {m : Nat}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes output : Array HenselNode} {target : SparsePolyZZ}
    (hcorrect : HenselLiftRecursiveCorrect termination m tree nodes target
      output) :
    HenselArrayReduces m nodes output := by
  induction hcorrect with
  | leaf index nodes stored target inputNode lifted parent hnode hstep
      hstepCorrect hstored hparent =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      apply henselArrayReduces_set m nodes index lifted hindex
      rw [hinput]
      exact ⟨hstepCorrect.2.1, hstepCorrect.2.2.1,
        hstepCorrect.2.2.2.1, hstepCorrect.2.2.2.2.1⟩
  | left index left nodes stored nodesAfterLeft target inputNode lifted parent
      hnode hstep hstepCorrect hstored hleftRun hleftCorrect hparent ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      exact (henselArrayReduces_set m nodes index lifted hindex (by
        rw [hinput]
        exact ⟨hstepCorrect.2.1, hstepCorrect.2.2.1,
          hstepCorrect.2.2.2.1, hstepCorrect.2.2.2.2.1⟩)).trans ih
  | right index right nodes stored output target inputNode lifted parent hnode
      hstep hstepCorrect hstored hparent hrightRun hrightCorrect ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      exact (henselArrayReduces_set m nodes index lifted hindex (by
        rw [hinput]
        exact ⟨hstepCorrect.2.1, hstepCorrect.2.2.1,
          hstepCorrect.2.2.2.1, hstepCorrect.2.2.2.2.1⟩)).trans ih
  | branch index left right nodes stored nodesAfterLeft output target inputNode
      lifted parent hnode hstep hstepCorrect hstored hleftRun hleftCorrect
      hparent hrightRun hrightCorrect leftIH rightIH =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index]! = inputNode := by
        have heq : nodes[index] = inputNode :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      exact (((henselArrayReduces_set m nodes index lifted hindex (by
        rw [hinput]
        exact ⟨hstepCorrect.2.1, hstepCorrect.2.2.1,
          hstepCorrect.2.2.2.1, hstepCorrect.2.2.2.2.1⟩)).trans leftIH).trans
            rightIH)

/-- The actual well-founded tree traversal preserves canonical sparse factor
representations at every array slot.  Each selected slot is replaced by the
canonical `g`/`h` returned by its concrete Hensel step; recursive calls then
apply the same frame argument to the updated array. -/
theorem HenselLiftRecursiveCorrect.arrayCanonical
    {termination : Generated.StrictHensel.DivmodTermination} {m : Nat}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes output : Array HenselNode} {target : SparsePolyZZ}
    (hcorrect : HenselLiftRecursiveCorrect termination m tree nodes target
      output)
    (hnodes : HenselArrayCanonical nodes) :
    HenselArrayCanonical output := by
  induction hcorrect with
  | leaf index nodes stored target inputNode lifted parent hnode hstep
      hstepCorrect hstored hparent =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      exact henselArrayCanonical_set nodes index lifted hindex hnodes
        ⟨hstepCorrect.2.2.2.2.2.1, hstepCorrect.2.2.2.2.2.2.1⟩
  | left index left nodes stored nodesAfterLeft target inputNode lifted parent
      hnode hstep hstepCorrect hstored hleftRun hleftCorrect hparent ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      exact ih (henselArrayCanonical_set nodes index lifted hindex hnodes
        ⟨hstepCorrect.2.2.2.2.2.1, hstepCorrect.2.2.2.2.2.2.1⟩)
  | right index right nodes stored output target inputNode lifted parent hnode
      hstep hstepCorrect hstored hparent hrightRun hrightCorrect ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      exact ih (henselArrayCanonical_set nodes index lifted hindex hnodes
        ⟨hstepCorrect.2.2.2.2.2.1, hstepCorrect.2.2.2.2.2.2.1⟩)
  | branch index left right nodes stored nodesAfterLeft output target inputNode
      lifted parent hnode hstep hstepCorrect hstored hleftRun hleftCorrect
      hparent hrightRun hrightCorrect leftIH rightIH =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      exact rightIH (leftIH (henselArrayCanonical_set nodes index lifted
        hindex hnodes
          ⟨hstepCorrect.2.2.2.2.2.1,
            hstepCorrect.2.2.2.2.2.2.1⟩))

/-- One complete execution of the generated recursive tree preserves a
physical coefficient-one head in every node's right factor.  Each updated
slot uses the concrete single-step head theorem; all other slots are framed
by the exact array write. -/
theorem HenselLiftRecursiveCorrect.arrayHOneHead
    {termination : Generated.StrictHensel.DivmodTermination} {m : Nat}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes output : Array HenselNode} {target : SparsePolyZZ}
    (hcorrect : HenselLiftRecursiveCorrect termination m tree nodes target
      output) (hm : 1 < m) (hnodes : HenselArrayHOneHead nodes) :
    HenselArrayHOneHead output := by
  induction hcorrect with
  | leaf index nodes stored target inputNode lifted parent hnode hstep
      hstepCorrect hstored hparent =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index] = inputNode :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      rcases hnodes index hindex with ⟨head, tail, hhead⟩
      have hlifted : HasPhysicalOneHead lifted.h := by
        rcases hstepCorrect.2.2.2.2.2.2.2 head tail (by
          simpa [hinput] using hhead) hm with ⟨suffix, hsuffix⟩
        exact ⟨head, suffix, hsuffix⟩
      exact henselArrayHOneHead_set nodes index lifted hindex hnodes hlifted
  | left index left nodes stored nodesAfterLeft target inputNode lifted parent
      hnode hstep hstepCorrect hstored hleftRun hleftCorrect hparent ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index] = inputNode :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      rcases hnodes index hindex with ⟨head, tail, hhead⟩
      have hlifted : HasPhysicalOneHead lifted.h := by
        rcases hstepCorrect.2.2.2.2.2.2.2 head tail (by
          simpa [hinput] using hhead) hm with ⟨suffix, hsuffix⟩
        exact ⟨head, suffix, hsuffix⟩
      exact ih
        (henselArrayHOneHead_set nodes index lifted hindex hnodes hlifted)
  | right index right nodes stored output target inputNode lifted parent hnode
      hstep hstepCorrect hstored hparent hrightRun hrightCorrect ih =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index] = inputNode :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      rcases hnodes index hindex with ⟨head, tail, hhead⟩
      have hlifted : HasPhysicalOneHead lifted.h := by
        rcases hstepCorrect.2.2.2.2.2.2.2 head tail (by
          simpa [hinput] using hhead) hm with ⟨suffix, hsuffix⟩
        exact ⟨head, suffix, hsuffix⟩
      exact ih
        (henselArrayHOneHead_set nodes index lifted hindex hnodes hlifted)
  | branch index left right nodes stored nodesAfterLeft output target inputNode
      lifted parent hnode hstep hstepCorrect hstored hleftRun hleftCorrect
      hparent hrightRun hrightCorrect leftIH rightIH =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hinput : nodes[index] = inputNode :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      rcases hnodes index hindex with ⟨head, tail, hhead⟩
      have hlifted : HasPhysicalOneHead lifted.h := by
        rcases hstepCorrect.2.2.2.2.2.2.2 head tail (by
          simpa [hinput] using hhead) hm with ⟨suffix, hsuffix⟩
        exact ⟨head, suffix, hsuffix⟩
      exact rightIH (leftIH
        (henselArrayHOneHead_set nodes index lifted hindex hnodes hlifted))

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

theorem HenselLiftLoopCorrect.topologyEq
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM) :
    HenselArrayTopologyEq initialNodes outputNodes := by
  induction hcorrect with
  | done => exact .refl _
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      exact hiteration.topologyEq.trans ih

/-- All quadratic rounds preserve every node modulo any divisor of the
initial modulus.  Instantiating `base = initialM = p` is the exact bridge from
the extracted lifted leaves back to the source finite-field factors. -/
theorem HenselLiftLoopCorrect.arrayReduces_of_dvd
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM base : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM)
    (hdiv : base ∣ initialM) :
    HenselArrayReduces base initialNodes outputNodes := by
  induction hcorrect generalizing base with
  | done m nodes hstop => exact HenselArrayReduces.refl base nodes
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      have hhead := hiteration.arrayReduces.of_dvd hdiv
      have hnextDiv : base ∣ m * m := dvd_mul_of_dvd_left hdiv m
      exact hhead.trans (ih hnextDiv)

/-- Every concrete quadratic lifting round preserves canonical factor fields
through the actual recursive tree execution recorded at that round. -/
theorem HenselLiftLoopCorrect.arrayCanonical
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM)
    (hinitial : HenselArrayCanonical initialNodes) :
    HenselArrayCanonical outputNodes := by
  induction hcorrect with
  | done => exact hinitial
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      exact ih (hiteration.arrayCanonical hinitial)

/-- Every actual quadratic round preserves physical monicity of all stored
right factors.  The lower bound on the initial modulus is propagated through
the source update `m := m * m`. -/
theorem HenselLiftLoopCorrect.arrayHOneHead
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM)
    (hm : 2 ≤ initialM) (hinitial : HenselArrayHOneHead initialNodes) :
    HenselArrayHOneHead outputNodes := by
  induction hcorrect with
  | done => exact hinitial
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      have hmOne : 1 < m := by omega
      have hnext : 2 ≤ m * m := by nlinarith
      exact ih hnext (hiteration.arrayHOneHead hmOne hinitial)

/-- The concrete quadratic loop only replaces its modulus by its square, so
the initial modulus divides the modulus returned by every actual execution. -/
theorem HenselLiftLoopCorrect.initialM_dvd_outputM
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM) :
    initialM ∣ outputM := by
  induction hcorrect with
  | done => exact dvd_refl _
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      exact dvd_trans ⟨m, by simp⟩ ih

/-- The modulus returned by the concrete quadratic Hensel loop is strictly
larger than the source target.  This is the precision fact consumed by
recombination; it follows from the actual stopping test rather than from an
independently supplied bound. -/
theorem HenselLiftLoopCorrect.outputM_gt_target
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM) :
    target < outputM := by
  induction hcorrect with
  | done m nodes hstop => omega
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      exact ih

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

/-- Extraction safety depends only on concrete node lookup results.  This
allows a certificate produced while constructing one subtree to be transported
across later array growth and disjoint node updates. -/
theorem HenselExtractInvariant.of_getElem?_eq
    (tree : Generated.StrictHensel.HenselLiftTree)
    (before after : Array HenselNode)
    (hlookup : ∀ index : Nat, after[index]? = before[index]?)
    (hinvariant : HenselExtractInvariant tree before) :
    HenselExtractInvariant tree after := by
  have harray : after = before := by
    apply Array.toList_inj.mp
    apply List.ext_getElem?
    intro index
    simpa only [Array.getElem?_toList] using hlookup index
  simpa [harray] using hinvariant

/-- An explicit, constructor-shaped extraction certificate.  Unlike the
recursor-normal form above, this relation supports ordinary induction and is
therefore convenient for proving preservation across array-prefix growth. -/
inductive HenselExtractCertificate (lower : Nat) :
    Generated.StrictHensel.HenselLiftTree → Array HenselNode → Prop
  | node
      (index : Nat)
      (left right : Option Generated.StrictHensel.HenselLiftTree)
      (nodes : Array HenselNode) (value : HenselNode)
      (hlower : lower ≤ index)
      (hnode : nodes[index]? = some value)
      (hleft : liftChildMatches value.left left)
      (hright : liftChildMatches value.right right)
      (hleftCertificate : ∀ child, left = some child →
        HenselExtractCertificate lower child nodes)
      (hrightCertificate : ∀ child, right = some child →
        HenselExtractCertificate lower child nodes) :
      HenselExtractCertificate lower (.node index left right) nodes

theorem HenselExtractCertificate.toInvariant
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes : Array HenselNode}
    {lower : Nat} (hcertificate : HenselExtractCertificate lower tree nodes) :
    HenselExtractInvariant tree nodes := by
  induction hcertificate with
  | node index left right nodes value hlower hnode hleft hright
      hleftCertificate hrightCertificate leftIH rightIH =>
      simp only [HenselExtractInvariant]
      refine ⟨value, hnode, hleft, hright, ?_, ?_⟩
      · cases left with
        | none => trivial
        | some child => exact leftIH child rfl
      · cases right with
        | none => trivial
        | some child => exact rightIH child rfl

theorem HenselExtractInvariant.of_topologyEq
    {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (htopology : HenselArrayTopologyEq before after)
    (hinvariant : HenselExtractInvariant tree before) :
    HenselExtractInvariant tree after := by
  exact Generated.StrictHensel.HenselLiftTree.rec
    (motive_1 := fun tree => ∀ before after,
      HenselArrayTopologyEq before after → HenselExtractInvariant tree before →
        HenselExtractInvariant tree after)
    (motive_2 := fun child => ∀ before after,
      HenselArrayTopologyEq before after →
      (child.elim True fun tree => HenselExtractInvariant tree before) →
      (child.elim True fun tree => HenselExtractInvariant tree after))
    (fun index left right leftIH rightIH before after htopology hinvariant => by
      simp only [HenselExtractInvariant] at hinvariant ⊢
      rcases hinvariant with
        ⟨value, hvalue, hleft, hright, hleftInv, hrightInv⟩
      have hbefore : index < before.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hvalue
        contradiction
      have hafter : index < after.size := by simpa [← htopology.1]
      have hbeforeValue : before[index] = value :=
        Option.some.inj ((Array.getElem?_eq_getElem hbefore).symm.trans hvalue)
      have hnodeTopology := htopology.2 index hbefore hafter
      refine ⟨after[index], Array.getElem?_eq_getElem hafter, ?_, ?_, ?_, ?_⟩
      · simpa [hbeforeValue, ← hnodeTopology.1] using hleft
      · simpa [hbeforeValue, ← hnodeTopology.2.1] using hright
      · cases left with
        | none => trivial
        | some child => exact leftIH before after htopology hleftInv
      · cases right with
        | none => trivial
        | some child => exact rightIH before after htopology hrightInv)
    (fun _ _ _ _ => trivial)
    (fun child childIH before after htopology hinvariant =>
      childIH before after htopology hinvariant)
    tree before after htopology hinvariant

/-- A constructor-shaped extraction certificate survives every array update
that preserves the structural fields at each concrete index. -/
theorem HenselExtractCertificate.of_topologyEq
    {lower : Nat} {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (htopology : HenselArrayTopologyEq before after)
    (hcertificate : HenselExtractCertificate lower tree before) :
    HenselExtractCertificate lower tree after := by
  induction hcertificate with
  | node index left right nodes value hlower hnode hleft hright
      hleftCertificate hrightCertificate leftIH rightIH =>
      have hbefore : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hafter : index < after.size := by simpa [← htopology.1]
      have hbeforeValue : nodes[index] = value :=
        Option.some.inj ((Array.getElem?_eq_getElem hbefore).symm.trans hnode)
      have hnodeTopology := htopology.2 index hbefore hafter
      refine .node index left right after after[index] hlower
        (Array.getElem?_eq_getElem hafter) ?_ ?_ ?_ ?_
      · have hfield : after[index].left = value.left := by
          simpa [hbeforeValue] using hnodeTopology.1.symm
        simpa [hfield] using hleft
      · have hfield : after[index].right = value.right := by
          simpa [hbeforeValue] using hnodeTopology.2.1.symm
        simpa [hfield] using hright
      · intro child hchild
        exact leftIH child hchild htopology
      · intro child hchild
        exact rightIH child hchild htopology

/-- Prefix preservation relation used by the generated append-only tree
builder: every lookup that succeeded before has the identical result after. -/
def HenselTreePrefix (before after : Array HenselNode) : Prop :=
  before.size ≤ after.size ∧
  ∀ index (hindex : index < before.size), after[index]? = before[index]?

def HenselTreePreservesFrom (lower : Nat)
    (before after : Array HenselNode) : Prop :=
  before.size ≤ after.size ∧
  ∀ index, lower ≤ index → index < before.size →
    after[index]? = before[index]?

theorem HenselExtractCertificate.of_prefix
    {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (hprefix : HenselTreePrefix before after)
    {lower : Nat} (hcertificate : HenselExtractCertificate lower tree before) :
    HenselExtractCertificate lower tree after := by
  induction hcertificate with
  | node index left right nodes value hlower hnode hleft hright
      hleftCertificate hrightCertificate leftIH rightIH =>
      have hindex : index < nodes.size := by
        by_contra hnot
        have : nodes[index]? = none :=
          Array.getElem?_eq_none (by omega)
        rw [this] at hnode
        contradiction
      exact .node index left right after value hlower
        ((hprefix.2 index hindex).trans hnode) hleft hright
        (fun child hchild => leftIH child hchild hprefix)
        (fun child hchild => rightIH child hchild hprefix)

theorem HenselExtractCertificate.of_preservesFrom
    {lower : Nat} {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (hpreserves : HenselTreePreservesFrom lower before after)
    (hcertificate : HenselExtractCertificate lower tree before) :
    HenselExtractCertificate lower tree after := by
  induction hcertificate with
  | node index left right nodes value hlower hnode hleft hright
      hleftCertificate hrightCertificate leftIH rightIH =>
      have hindex : index < nodes.size := by
        by_contra hnot
        have : nodes[index]? = none := Array.getElem?_eq_none (by omega)
        rw [this] at hnode
        contradiction
      exact .node index left right after value hlower
        ((hpreserves.2 index hlower hindex).trans hnode) hleft hright
        (fun child hchild => leftIH child hchild hpreserves)
        (fun child hchild => rightIH child hchild hpreserves)

theorem HenselExtractCertificate.lower_mono
    {lower lower' : Nat} {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes : Array HenselNode} (hlower : lower' ≤ lower)
    (hcertificate : HenselExtractCertificate lower tree nodes) :
    HenselExtractCertificate lower' tree nodes := by
  induction hcertificate with
  | node index left right nodes value hindexLower hnode hleft hright
      hleftCertificate hrightCertificate leftIH rightIH =>
      exact .node index left right nodes value (by omega) hnode hleft hright
        (fun child hchild => leftIH child hchild)
        (fun child hchild => rightIH child hchild)

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

def HenselFactorArrayCanonical (factors : Array SparsePolyZZ) : Prop :=
  ∀ factor ∈ factors.toList,
    StrictPolynomialMod.SparsePolyZZCanonical factor

theorem HenselFactorArrayCanonical.push
    {factors : Array SparsePolyZZ} {factor : SparsePolyZZ}
    (hfactors : HenselFactorArrayCanonical factors)
    (hfactor : StrictPolynomialMod.SparsePolyZZCanonical factor) :
    HenselFactorArrayCanonical (factors.push factor) := by
  intro candidate hcandidate
  simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hcandidate
  rcases hcandidate with hcandidate | rfl
  · exact hfactors candidate hcandidate
  · exact hfactor

/-- The exact extraction traversal appends only `g`/`h` fields read from the
canonical lifted node array, so it preserves canonicality of the input prefix
and establishes it for every newly extracted factor. -/
theorem HenselExtractCorrect.outputCanonical
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes : Array HenselNode} {input output : Array SparsePolyZZ}
    (hcorrect : HenselExtractCorrect tree nodes input output)
    (hnodes : HenselArrayCanonical nodes)
    (hinput : HenselFactorArrayCanonical input) :
    HenselFactorArrayCanonical output := by
  induction hcorrect with
  | leaf index nodes input output node hnode houtput =>
      subst output
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlookup : nodes[index] = node :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      have hnodeCanonical : HenselNodeCanonical node := by
        simpa [hlookup] using hnodes index hindex
      exact (hinput.push hnodeCanonical.1).push hnodeCanonical.2
  | left index left nodes input afterLeft output node hnode hleftRun
      hleftCorrect houtput ih =>
      subst output
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlookup : nodes[index] = node :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      have hnodeCanonical : HenselNodeCanonical node := by
        simpa [hlookup] using hnodes index hindex
      exact (ih hnodes hinput).push hnodeCanonical.2
  | right index right nodes input afterLeft output node hnode hafterLeft
      hrightRun hrightCorrect ih =>
      subst afterLeft
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlookup : nodes[index] = node :=
        Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
      have hnodeCanonical : HenselNodeCanonical node := by
        simpa [hlookup] using hnodes index hindex
      exact ih hnodes (hinput.push hnodeCanonical.1)
  | branch index left right nodes input afterLeft output node hnode hleftRun
      hleftCorrect hrightRun hrightCorrect leftIH rightIH =>
      exact rightIH hnodes (leftIH hnodes hinput)

/-- Pure denotation of the exact leaf order used by the generated extraction
walk.  It reads the same node `g`/`h` fields selected by the source branches;
it does not construct or normalize any factor. -/
def henselExtractedFactors :
    Generated.StrictHensel.HenselLiftTree → Array HenselNode →
      List SparsePolyZZ
  | .node index left right, nodes =>
      (match left with
        | none => [nodes[index]!.g]
        | some child => henselExtractedFactors child nodes) ++
      (match right with
        | none => [nodes[index]!.h]
        | some child => henselExtractedFactors child nodes)
termination_by tree _ => tree.nodeCount
decreasing_by
  · exact Generated.StrictHensel.HenselLiftTree.left_nodeCount_lt
      index child right
  · exact Generated.StrictHensel.HenselLiftTree.right_nodeCount_lt
      index left child

/-- The pure extraction denotation depends only on the concrete array cells
visited by its topology.  This is the exact transport lemma used when the
source-order right traversal preserves all indices of the completed left
subtree. -/
theorem henselExtractedFactors_eq_of_lookups
    (tree : Generated.StrictHensel.HenselLiftTree)
    (before after : Array HenselNode)
    (hlookup : ∀ index ∈ henselLiftTreeIndices tree,
      after[index]? = before[index]?) :
    henselExtractedFactors tree after =
      henselExtractedFactors tree before := by
  cases htree : tree with
  | node index left right =>
      rw [htree] at hlookup
      have hroot : after[index]? = before[index]? :=
        hlookup index (by
          cases left <;> cases right <;>
            simp [henselLiftTreeIndices])
      have hrootBang : after[index]! = before[index]! := by
        rw [getElem!_def, getElem!_def, hroot]
      rw [henselExtractedFactors.eq_def, henselExtractedFactors.eq_def]
      cases left with
      | none =>
          cases right with
          | none => simp [hrootBang]
          | some child =>
              have hright := henselExtractedFactors_eq_of_lookups child
                before after (by
                  intro childIndex hchildIndex
                  exact hlookup childIndex (by
                    rw [henselLiftTreeIndices.eq_def]
                    simp [hchildIndex]))
              simp [hrootBang, hright]
      | some child =>
          cases right with
          | none =>
              have hleft := henselExtractedFactors_eq_of_lookups child
                before after (by
                  intro childIndex hchildIndex
                  exact hlookup childIndex (by
                    rw [henselLiftTreeIndices.eq_def]
                    simp [hchildIndex]))
              simp [hrootBang, hleft]
          | some rightChild =>
              have hleft := henselExtractedFactors_eq_of_lookups child
                before after (by
                  intro childIndex hchildIndex
                  exact hlookup childIndex (by
                    rw [henselLiftTreeIndices.eq_def]
                    simp [hchildIndex]))
              have hright := henselExtractedFactors_eq_of_lookups rightChild
                before after (by
                  intro childIndex hchildIndex
                  exact hlookup childIndex (by
                    rw [henselLiftTreeIndices.eq_def]
                    simp [hchildIndex]))
              simp [hleft, hright]
termination_by tree.nodeCount
decreasing_by
  all_goals subst tree
  all_goals simp_all [Generated.StrictHensel.HenselLiftTree.nodeCount]
  all_goals omega

/-- The product of the leaves left by one actual recursive Hensel traversal is
the traversal target modulo the exact squared source modulus.  The `Nodup`
premise is discharged for the generated canonical tree by
`henselTreeBuildTopology_indices_nodup_bounded`; it is used here to justify
the real mutable-array frames between the source-order child calls. -/
theorem HenselLiftRecursiveCorrect.extractedFactors_product
    {termination : Generated.StrictHensel.DivmodTermination} {m : Nat}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes output : Array HenselNode} {target : SparsePolyZZ}
    (hcorrect : HenselLiftRecursiveCorrect termination m tree nodes target
      output)
    (hnodup : (henselLiftTreeIndices tree).Nodup) :
    ((henselExtractedFactors tree output).map
      (toPolyMod (m ^ 2))).prod = toPolyMod (m ^ 2) target := by
  induction hcorrect with
  | leaf index nodes stored target inputNode lifted parent hnode hstep
      hstepCorrect hstored hparent =>
      subst stored
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlifted : (nodes.set! index lifted)[index]? = some lifted := by
        simp [Array.set!, hindex]
      have hparentEq : parent = lifted :=
        Option.some.inj (hparent.symm.trans hlifted)
      subst parent
      have hbang : (nodes.set! index lifted)[index]! = lifted := by
        simp [Array.set!, hindex]
      have hbang' : (nodes.setIfInBounds index lifted)[index]! = lifted := by
        simpa only [Array.set!_eq_setIfInBounds] using hbang
      simpa [henselExtractedFactors, hbang'] using hstepCorrect.1.1
  | left index left nodes stored nodesAfterLeft target inputNode lifted parent
      hnode hstep hstepCorrect hstored hleftRun hleftCorrect hparent ih =>
      subst stored
      have hparts : index ∉ henselLiftTreeIndices left ∧
          (henselLiftTreeIndices left).Nodup := by
        simpa [henselLiftTreeIndices] using hnodup
      have hrootPreserved := hleftCorrect.preserves_not_mem index hparts.1
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlifted : (nodes.set! index lifted)[index]? = some lifted := by
        simp [Array.set!, hindex]
      have hparentEq : parent = lifted :=
        Option.some.inj (hparent.symm.trans (hrootPreserved.trans hlifted))
      subst parent
      have hbang : nodesAfterLeft[index]! = lifted := by
        rw [getElem!_def, hparent]
      have hleftProduct := ih hparts.2
      simpa [henselExtractedFactors, hbang, hleftProduct] using
        hstepCorrect.1.1
  | right index right nodes stored output target inputNode lifted parent hnode
      hstep hstepCorrect hstored hparent hrightRun hrightCorrect ih =>
      subst stored
      have hparts : index ∉ henselLiftTreeIndices right ∧
          (henselLiftTreeIndices right).Nodup := by
        simpa [henselLiftTreeIndices] using hnodup
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlifted : (nodes.set! index lifted)[index]? = some lifted := by
        simp [Array.set!, hindex]
      have hparentEq : parent = lifted :=
        Option.some.inj (hparent.symm.trans hlifted)
      subst parent
      have hrootPreserved := hrightCorrect.preserves_not_mem index hparts.1
      have hbang : output[index]! = lifted := by
        rw [getElem!_def, hrootPreserved, hlifted]
      have hrightProduct := ih hparts.2
      simpa [henselExtractedFactors, hbang, hrightProduct] using
        hstepCorrect.1.1
  | branch index left right nodes stored nodesAfterLeft output target inputNode
      lifted parent hnode hstep hstepCorrect hstored hleftRun hleftCorrect
      hparent hrightRun hrightCorrect leftIH rightIH =>
      subst stored
      have hparts : index ∉ henselLiftTreeIndices left ∧
          index ∉ henselLiftTreeIndices right ∧
          (henselLiftTreeIndices left).Nodup ∧
          (henselLiftTreeIndices right).Nodup ∧
          List.Disjoint (henselLiftTreeIndices left)
            (henselLiftTreeIndices right) := by
        have hparts :
            (index ∉ henselLiftTreeIndices left ∧
              index ∉ henselLiftTreeIndices right) ∧
            (henselLiftTreeIndices left).Nodup ∧
            (henselLiftTreeIndices right).Nodup ∧
            ∀ a ∈ henselLiftTreeIndices left,
              ∀ b ∈ henselLiftTreeIndices right, ¬a = b := by
          simpa [henselLiftTreeIndices, List.nodup_append] using hnodup
        exact ⟨hparts.1.1, hparts.1.2, hparts.2.1, hparts.2.2.1,
          by simpa [List.disjoint_iff_ne] using hparts.2.2.2⟩
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlifted : (nodes.set! index lifted)[index]? = some lifted := by
        simp [Array.set!, hindex]
      have hrootAfterLeft := hleftCorrect.preserves_not_mem index hparts.1
      have hparentEq : parent = lifted :=
        Option.some.inj (hparent.symm.trans (hrootAfterLeft.trans hlifted))
      subst parent
      have hleftFinal : henselExtractedFactors left output =
          henselExtractedFactors left nodesAfterLeft := by
        apply henselExtractedFactors_eq_of_lookups
        intro childIndex hchild
        exact hrightCorrect.preserves_not_mem childIndex (by
          intro hright
          exact hparts.2.2.2.2 hchild hright)
      have hleftProduct := leftIH hparts.2.2.1
      have hrightProduct := rightIH hparts.2.2.2.1
      simpa [henselExtractedFactors, hleftFinal, hleftProduct,
        hrightProduct] using hstepCorrect.1.1

/-- The actual quadratic-precision loop carries the concrete leaf-product
invariant from its initial modulus to the exact modulus returned by the
generated well-founded loop.  Each recursive edge is justified by the
corresponding generated tree traversal at the current modulus. -/
theorem HenselLiftLoopCorrect.extractedFactors_product
    {termination : Generated.StrictHensel.DivmodTermination}
    {tree : Generated.StrictHensel.HenselLiftTree} {f : SparsePolyZZ}
    {target initialM outputM : Nat}
    {initialNodes outputNodes : Array HenselNode}
    (hcorrect : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes outputNodes outputM)
    (hnodup : (henselLiftTreeIndices tree).Nodup)
    (hinitial : ((henselExtractedFactors tree initialNodes).map
      (toPolyMod initialM)).prod = toPolyMod initialM f) :
    ((henselExtractedFactors tree outputNodes).map
      (toPolyMod outputM)).prod = toPolyMod outputM f := by
  induction hcorrect with
  | done => exact hinitial
  | step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
      htail ih =>
      apply ih
      rw [← Nat.pow_two]
      exact hiteration.extractedFactors_product hnodup

/-- The semantic extraction trace returns the input prefix followed by exactly
the node fields named by `henselExtractedFactors`, in source left-to-right
order. -/
theorem HenselExtractCorrect.toList_eq
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes : Array HenselNode} {input output : Array SparsePolyZZ}
    (hcorrect : HenselExtractCorrect tree nodes input output) :
    output.toList = input.toList ++ henselExtractedFactors tree nodes := by
  induction hcorrect with
  | leaf index nodes input output node hnode houtput =>
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlookup : nodes[index]! = node := by
        have heq : nodes[index] = node :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      subst output
      simp [henselExtractedFactors, hlookup, List.append_assoc]
  | left index left nodes input afterLeft output node hnode hleftRun
      hleftCorrect houtput ih =>
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlookup : nodes[index]! = node := by
        have heq : nodes[index] = node :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      subst output
      simp [henselExtractedFactors, hlookup, ih, List.append_assoc]
  | right index right nodes input afterLeft output node hnode hafterLeft
      hrightRun hrightCorrect ih =>
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hlookup : nodes[index]! = node := by
        have heq : nodes[index] = node :=
          Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
      subst afterLeft
      simpa [henselExtractedFactors, hlookup, List.append_assoc] using ih
  | branch index left right nodes input afterLeft output node hnode hleftRun
      hleftCorrect hrightRun hrightCorrect leftIH rightIH =>
      rw [rightIH, leftIH]
      simp [henselExtractedFactors, List.append_assoc]

/-- Reading the same extraction topology from two node arrays related by
`HenselArrayReduces` produces pointwise congruent factor lists. -/
theorem HenselExtractCertificate.extractedFactors_forall₂
    {lower p : Nat} {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (hcertificate : HenselExtractCertificate lower tree before)
    (hreduce : HenselArrayReduces p before after) :
    List.Forall₂ (fun old new => toPolyMod p new = toPolyMod p old)
      (henselExtractedFactors tree before)
      (henselExtractedFactors tree after) := by
  induction hcertificate with
  | node index left right nodes value hlower hnode hleft hright
      hleftCertificate hrightCertificate leftIH rightIH =>
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hafterIndex : index < after.size := by
        rw [← hreduce.1]
        exact hindex
      have hnodeReduce := hreduce.2 index hindex hafterIndex
      have hnodeReduceBang :
          HenselNodeReduces p nodes[index]! after[index]! := by
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex,
          Array.getElem?_eq_getElem hafterIndex] using hnodeReduce
      rw [henselExtractedFactors.eq_def, henselExtractedFactors.eq_def]
      apply List.rel_append
      · cases left with
        | none =>
            exact .cons hnodeReduceBang.1 .nil
        | some child =>
            exact leftIH child rfl hreduce
      · cases right with
        | none =>
            exact .cons hnodeReduceBang.2.1 .nil
        | some child =>
            exact rightIH child rfl hreduce

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

/-- Genuine composition of the actual well-founded lift-loop result with the
actual generated extraction traversal.  The builder certificate is transported
through the concrete loop trace using structural-field preservation; no tree
or extracted factor array is supplied by a specification oracle. -/
theorem henselLiftLoopThenExtractRawIR_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes liftedNodes : Array HenselNode)
    (outputM : Nat) (factors : Array SparsePolyZZ)
    (hloop : HenselLiftLoopCorrect termination tree f target initialM
      initialNodes liftedNodes outputM)
    (hcertificate : HenselExtractCertificate 0 tree initialNodes) :
    ∃ output,
      Generated.StrictHensel.__hensel_extract_factors_raw_ir tree liftedNodes
          factors = .ok output ∧
      HenselExtractCorrect tree liftedNodes factors output := by
  have hliftedCertificate :=
    hcertificate.of_topologyEq hloop.topologyEq
  exact __hensel_extract_factors_raw_ir_refines tree liftedNodes factors
    hliftedCertificate.toInvariant

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

/-- Pointwise semantic effect of the final source normalization: each output
factor differs from the corresponding input by a unit scalar modulo `p`. -/
def HenselNormalizeUnitRel (p : Nat)
    (before after : Array SparsePolyZZ) : Prop :=
  before.size = after.size ∧
  ∀ index (hbefore : index < before.size) (hafter : index < after.size),
    ∃ scale : ZMod p, IsUnit scale ∧
      toPolyMod p after[index] =
        Polynomial.C scale * toPolyMod p before[index]

theorem HenselNormalizeUnitRel.irreducible
    {p : Nat} {before after : Array SparsePolyZZ}
    (hrel : HenselNormalizeUnitRel p before after)
    (hirreducible : ∀ index (hindex : index < before.size),
      Irreducible (toPolyMod p before[index])) :
    ∀ index (hindex : index < after.size),
      Irreducible (toPolyMod p after[index]) := by
  intro index hafter
  have hbefore : index < before.size := by rw [hrel.1]; exact hafter
  rcases hrel.2 index hbefore hafter with ⟨scale, hscale, heq⟩
  rw [heq]
  exact (irreducible_isUnit_mul (Polynomial.isUnit_C.mpr hscale)).2
    (hirreducible index hbefore)

private theorem isCoprime_unit_mul_unit_mul {R : Type*} [CommRing R]
    {left right leftUnit rightUnit : R}
    (hleftUnit : IsUnit leftUnit) (hrightUnit : IsUnit rightUnit)
    (hcoprime : IsCoprime left right) :
    IsCoprime (leftUnit * left) (rightUnit * right) := by
  rcases hleftUnit with ⟨leftUnit, rfl⟩
  rcases hrightUnit with ⟨rightUnit, rfl⟩
  rcases hcoprime with ⟨a, b, hab⟩
  refine ⟨(↑leftUnit⁻¹ : R) * a, (↑rightUnit⁻¹ : R) * b, ?_⟩
  calc
    (↑leftUnit⁻¹ : R) * a * (↑leftUnit * left) +
        (↑rightUnit⁻¹ : R) * b * (↑rightUnit * right) =
      ((↑leftUnit⁻¹ : R) * ↑leftUnit) * (a * left) +
        ((↑rightUnit⁻¹ : R) * ↑rightUnit) * (b * right) := by ring
    _ = a * left + b * right := by simp
    _ = 1 := hab

/-- Pairwise coprimality is invariant under the exact pointwise unit scaling
performed by the generated final Hensel normalization block. -/
theorem HenselNormalizeUnitRel.isCoprime
    {p : Nat} {before after : Array SparsePolyZZ}
    (hrel : HenselNormalizeUnitRel p before after)
    {left right : Nat} (hleft : left < after.size)
    (hright : right < after.size)
    (hcoprime : IsCoprime (toPolyMod p before[left]!)
      (toPolyMod p before[right]!)) :
    IsCoprime (toPolyMod p after[left]!) (toPolyMod p after[right]!) := by
  have hleftBefore : left < before.size := by rw [hrel.1]; exact hleft
  have hrightBefore : right < before.size := by rw [hrel.1]; exact hright
  rcases hrel.2 left hleftBefore hleft with
    ⟨leftScale, hleftScale, hleftEq⟩
  rcases hrel.2 right hrightBefore hright with
    ⟨rightScale, hrightScale, hrightEq⟩
  rw [getElem!_pos after left hleft,
    getElem!_pos after right hright]
  rw [hleftEq, hrightEq]
  have hcoprime' : IsCoprime (toPolyMod p before[left])
      (toPolyMod p before[right]) := by
    simpa [getElem!_pos before left hleftBefore,
      getElem!_pos before right hrightBefore] using hcoprime
  exact isCoprime_unit_mul_unit_mul
    (Polynomial.isUnit_C.mpr hleftScale)
    (Polynomial.isUnit_C.mpr hrightScale) hcoprime'

theorem HenselNormalizeCorrect.unitRel
    {p m : Nat} {before after : Array SparsePolyZZ}
    (hdiv : p ∣ m)
    (hcorrect : HenselNormalizeCorrect before (m : Int) after) :
    HenselNormalizeUnitRel p before after := by
  cases hcorrect with
  | empty hresult =>
      exact ⟨rfl, fun index hbefore hafter =>
        ⟨1, isUnit_one, by simp⟩⟩
  | alreadyOne first leading hresult hfirst hone =>
      exact ⟨rfl, fun index hbefore hafter =>
        ⟨1, isUnit_one, by simp⟩⟩
  | normalized first leading inverse scaled normalized hresult hfirst hnotOne
      hinverse hscaled hnormalized =>
      have hzero : 0 < before.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hresult
        contradiction
      have hfirstGet : before[0] = first :=
        Option.some.inj ((Array.getElem?_eq_getElem hzero).symm.trans hresult)
      have hinverseSuccess :
          (ZZ.invert 0 leading.2 (m : Int)).1 = true := by
        rw [hinverse]
      have hinverseM :
          (leading.2 : ZMod m) * (inverse : ZMod m) = 1 := by
        simpa [hinverse] using
          (CLPoly.Math.ZZ.invert_success_mul_eq_one leading.2 m
            hinverseSuccess)
      have hinverseP : (leading.2 : ZMod p) * (inverse : ZMod p) = 1 := by
        have hprojected := congrArg
          (ZMod.castHom hdiv (ZMod p)) hinverseM
        simpa only [map_mul, map_one, ZMod.castHom_apply,
          ZMod.cast_intCast hdiv] using hprojected
      have hinverseUnit : IsUnit (inverse : ZMod p) :=
        isUnit_iff_exists_inv.mpr ⟨(leading.2 : ZMod p), by
          simpa [mul_comm] using hinverseP⟩
      have hnormalizedSemantic :=
        __upoly_mod_coeff_raw_ir_preserves_divisor_of_run hdiv scaled
          normalized hnormalized
      rw [hscaled, scaleCoeffs_toPolyMod] at hnormalizedSemantic
      refine ⟨by simp [Array.set!, hzero], ?_⟩
      intro index hbefore hafter
      by_cases hindex : index = 0
      · subst index
        refine ⟨(inverse : ZMod p), hinverseUnit, ?_⟩
        simpa [Array.set!, hzero, hfirstGet] using hnormalizedSemantic
      · refine ⟨1, isUnit_one, ?_⟩
        change toPolyMod p (before.setIfInBounds 0 normalized)[index] = _
        have hzeroIndex : 0 ≠ index := by omega
        have hget : (before.setIfInBounds 0 normalized)[index] =
            before[index] := by
          simpa [hzeroIndex] using
            (@Array.getElem_setIfInBounds _ before 0 normalized index hbefore)
        rw [hget]
        simp

/-- Product-level form of the exact generated normalization trace.  The
source block either leaves the array unchanged or scales only its first
factor by the concrete modular inverse that it computed. -/
theorem HenselNormalizeCorrect.product_eq_unit_mul
    {m : Nat} {before after : Array SparsePolyZZ}
    (hcorrect : HenselNormalizeCorrect before (m : Int) after) :
    ∃ scale : ZMod m, IsUnit scale ∧
      (after.toList.map (toPolyMod m)).prod =
        Polynomial.C scale * (before.toList.map (toPolyMod m)).prod := by
  cases hcorrect with
  | empty hresult =>
      exact ⟨1, isUnit_one, by simp⟩
  | alreadyOne first leading hresult hfirst hone =>
      exact ⟨1, isUnit_one, by simp⟩
  | normalized first leading inverse scaled normalized hresult hfirst hnotOne
      hinverse hscaled hnormalized =>
      have hzero : 0 < before.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hresult
        contradiction
      have hfirstGet : before[0] = first :=
        Option.some.inj ((Array.getElem?_eq_getElem hzero).symm.trans hresult)
      have hinverseSuccess :
          (ZZ.invert 0 leading.2 (m : Int)).1 = true := by
        rw [hinverse]
      have hinverseM :
          (leading.2 : ZMod m) * (inverse : ZMod m) = 1 := by
        simpa [hinverse] using
          (CLPoly.Math.ZZ.invert_success_mul_eq_one leading.2 m
            hinverseSuccess)
      have hinverseUnit : IsUnit (inverse : ZMod m) :=
        isUnit_iff_exists_inv.mpr ⟨(leading.2 : ZMod m), by
          simpa [mul_comm] using hinverseM⟩
      have hnormalizedSemantic :=
        __upoly_mod_coeff_raw_ir_preserves_divisor_of_run (dvd_refl m)
          scaled normalized hnormalized
      rw [hscaled, scaleCoeffs_toPolyMod] at hnormalizedSemantic
      refine ⟨(inverse : ZMod m), hinverseUnit, ?_⟩
      change (List.map (toPolyMod m)
        (before.setIfInBounds 0 normalized).toList).prod = _
      rw [Array.toList_setIfInBounds]
      have hlist : before.toList = first :: before.toList.drop 1 := by
        rw [← hfirstGet]
        simpa using (List.getElem_cons_drop (as := before.toList) (i := 0)
          (by simpa using hzero)).symm
      rw [hlist]
      simp only [List.set_cons_zero, List.map_cons, List.prod_cons]
      rw [hnormalizedSemantic]
      ring

/-- The exact generated final-normalization branch preserves the sparse
integer representation invariant pointwise.  Only index zero can be
rewritten, and its replacement is the concrete scale/filter output. -/
theorem HenselNormalizeCorrect.canonical
    {before after : Array SparsePolyZZ} {m : ZZ}
    (hcorrect : HenselNormalizeCorrect before m after)
    (hbefore : ∀ index (hindex : index < before.size),
      StrictPolynomialMod.SparsePolyZZCanonical before[index]) :
    ∀ index (hindex : index < after.size),
      StrictPolynomialMod.SparsePolyZZCanonical after[index] := by
  cases hcorrect with
  | empty hresult => exact hbefore
  | alreadyOne first leading hresult hfirst hone => exact hbefore
  | normalized first leading inverse scaled normalized hresult hfirst hnotOne
      hinverse hscaled hnormalized =>
      intro index hafter
      have hzero : 0 < before.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hresult
        contradiction
      have hbeforeIndex : index < before.size := by simpa using hafter
      by_cases hindex : index = 0
      · subst index
        have hfirstGet : before[0] = first :=
          Option.some.inj
            ((Array.getElem?_eq_getElem hzero).symm.trans hresult)
        have hfirstCanonical := hbefore 0 hzero
        rw [hfirstGet] at hfirstCanonical
        have hnormalizedCanonical :=
          modCoeffOutput_scaleCoeffs_canonical first inverse m hfirstCanonical
        have hnormalizedEq : normalized =
            Generated.StrictHensel.modCoeffOutput
              (Generated.StrictHensel.scaleCoeffs first inverse) m := by
          rw [hscaled] at hnormalized
          simpa [Generated.StrictHensel.__upoly_mod_coeff_raw_ir] using
            (Except.ok.inj hnormalized).symm
        simpa [Array.set!, hzero, hnormalizedEq] using hnormalizedCanonical
      · change StrictPolynomialMod.SparsePolyZZCanonical
          (before.setIfInBounds 0 normalized)[index]
        have hget : (before.setIfInBounds 0 normalized)[index] =
            before[index] := by
          simpa [Ne.symm hindex] using
            (@Array.getElem_setIfInBounds _ before 0 normalized index
              hbeforeIndex)
        rw [hget]
        exact hbefore index hbeforeIndex

/-- Array-level form of `HenselNormalizeCorrect.canonical`, matching the
canonicality contract carried by extraction and consumed by recombination. -/
theorem HenselNormalizeCorrect.outputCanonical
    {before after : Array SparsePolyZZ} {m : ZZ}
    (hcorrect : HenselNormalizeCorrect before m after)
    (hbefore : HenselFactorArrayCanonical before) :
    HenselFactorArrayCanonical after := by
  intro factor hfactor
  rcases List.mem_iff_get.mp hfactor with ⟨index, hget⟩
  have hindex : index.1 < after.size := by simpa using index.2
  have hget' : after[index.1] = factor := by simpa using hget
  rw [← hget']
  exact hcorrect.canonical
    (fun beforeIndex hbeforeIndex =>
      hbefore before[beforeIndex] (Array.getElem_mem_toList hbeforeIndex))
    index.1 hindex

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

/-- The generated leading-coefficient baking write changes only array slot
zero, so every positive-index factor is definitionally the original one. -/
theorem HenselAdjustFirstFactorCorrect.getElem_eq_of_pos
    {f : SparsePolyZZ} {factors adjusted : Array SparsePolyZp} {p : UInt64}
    (hcorrect : HenselAdjustFirstFactorCorrect f factors p adjusted)
    (index : Nat) (hindex : index < adjusted.size) (hpos : 0 < index) :
    ∃ hsource : index < factors.size,
      getElem adjusted index hindex = getElem factors index hsource := by
  cases hcorrect with
  | adjusted leading first value hsource hfirst hvalue =>
      have hzero : 0 < factors.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hfirst
        contradiction
      have hsourceIndex : index < factors.size := by
        simpa using hindex
      refine ⟨hsourceIndex, ?_⟩
      simpa only [Array.set!_eq_setIfInBounds] using
        (Array.getElem_setIfInBounds_ne (xs := factors) (i := 0)
          (a := value) (j := index) hsourceIndex (by omega))

/-- Hence monicity of the finite-field factors produced by `__factor_Zp` is
preserved at every positive slot of the concrete adjusted array. -/
theorem HenselAdjustFirstFactorCorrect.monic_of_pos
    {f : SparsePolyZZ} {factors adjusted : Array SparsePolyZp} {p : UInt64}
    (hcorrect : HenselAdjustFirstFactorCorrect f factors p adjusted)
    (hmonic : ∀ index (hindex : index < factors.size),
      (CLPoly.Math.SparsePolyZp.toPoly p.toNat
        (getElem factors index hindex)).Monic)
    (index : Nat) (hindex : index < adjusted.size) (hpos : 0 < index) :
    (CLPoly.Math.SparsePolyZp.toPoly p.toNat
      (getElem adjusted index hindex)).Monic := by
  rcases hcorrect.getElem_eq_of_pos index hindex hpos with
    ⟨hsource, heq⟩
  rw [heq]
  exact hmonic index hsource

private theorem henselAdjustNormalization_toPoly
    (p : Nat) (factor : SparsePolyZp) :
    CLPoly.Math.SparsePolyZp.toPoly p (SparsePolyZp.normalization factor) =
      CLPoly.Math.SparsePolyZp.toPoly p factor := by
  unfold SparsePolyZp.normalization CLPoly.Math.SparsePolyZp.toPoly
  rw [Array.toList_filter]
  induction factor.toList with
  | nil => simp [CLPoly.Math.listSum]
  | cons term rest ih =>
      rcases term with ⟨monomial, coefficient⟩
      by_cases hzero : coefficient.val = 0
      · simpa [hzero, CLPoly.Math.listSum, CLPoly.Math.Zp.toZMod] using ih
      · simp [hzero, CLPoly.Math.listSum, ih]

private theorem scaleZpCoeffsList_toPoly
    (p : Nat) (hp2 : p * p ≤ UInt64.size) (coefficient : Zp)
    (hcoefficient : CLPoly.Math.Zp.Reduced p coefficient) :
    ∀ terms : List (UMonomial × Zp),
      CLPoly.Math.SparsePolyZp.AllReduced p terms →
      CLPoly.Math.listSum p
          (terms.map fun term => (term.1, term.2 * coefficient)) =
        Polynomial.C (CLPoly.Math.Zp.toZMod p coefficient) *
          CLPoly.Math.listSum p terms := by
  intro terms hterms
  induction terms with
  | nil => simp [CLPoly.Math.listSum]
  | cons term rest ih =>
      have hterm := hterms term List.mem_cons_self
      have hrest : CLPoly.Math.SparsePolyZp.AllReduced p rest :=
        fun item hitem => hterms item (List.mem_cons_of_mem term hitem)
      have hp : 0 < p := Nat.zero_lt_of_lt hterm.2
      have hnoOverflow :
          term.2.val.toNat * coefficient.val.toNat < UInt64.size := by
        calc
          term.2.val.toNat * coefficient.val.toNat < p * p :=
            Nat.mul_lt_mul_of_lt_of_le hterm.2
              (Nat.le_of_lt hcoefficient.2) hp
          _ ≤ UInt64.size := hp2
      have hmul := CLPoly.Math.Zp.toZMod_mul p term.2 coefficient
        hterm.1 hcoefficient.1 hnoOverflow
      rcases term with ⟨monomial, value⟩
      simp only [List.map_cons, CLPoly.Math.listSum_cons]
      rw [hmul, ih hrest, mul_add, Polynomial.C_mul_monomial]
      ring

/-- The exact coefficient map used by the source first-factor adjustment is
multiplication by the represented field constant. -/
theorem scaleZpCoeffs_toPoly
    (p : Nat) (hp2 : p * p ≤ UInt64.size) (coefficient : Zp)
    (hcoefficient : CLPoly.Math.Zp.Reduced p coefficient)
    (factor : SparsePolyZp)
    (hfactor : CLPoly.Math.SparsePolyZp.AllReduced p factor.toList) :
    CLPoly.Math.SparsePolyZp.toPoly p
        (Generated.StrictHensel.scaleZpCoeffs factor coefficient) =
      Polynomial.C (CLPoly.Math.Zp.toZMod p coefficient) *
        CLPoly.Math.SparsePolyZp.toPoly p factor := by
  unfold Generated.StrictHensel.scaleZpCoeffs
    CLPoly.Math.SparsePolyZp.toPoly
  rw [Array.toList_map]
  exact scaleZpCoeffsList_toPoly p hp2 coefficient hcoefficient
    factor.toList hfactor

theorem zpOfInt_reduced (coefficient : Int) (p : UInt64)
    (hp : 0 < p.toNat) :
    CLPoly.Math.Zp.Reduced p.toNat (Zp.ofInt coefficient p) := by
  constructor
  · simp [Zp.ofInt]
  · unfold Zp.ofInt
    have hpInt : (0 : Int) < p.toNat := by exact_mod_cast hp
    have hnonneg : 0 ≤ coefficient.emod p.toNat :=
      Int.emod_nonneg _ (by omega)
    have hltInt : coefficient.emod p.toNat < p.toNat :=
      Int.emod_lt_of_pos _ hpInt
    have hltNatP : (coefficient.emod p.toNat).toNat < p.toNat :=
      (Int.toNat_lt hnonneg).2 hltInt
    have hltSize : (coefficient.emod p.toNat).toNat < UInt64.size :=
      lt_trans hltNatP (UInt64.toNat_lt_size p)
    have hword : ((coefficient.emod p.toNat).toNat.toUInt64).toNat =
        (coefficient.emod p.toNat).toNat :=
      UInt64.toNat_ofNat_of_lt hltSize
    have hremNonneg : ¬coefficient % (p.toNat : Int) < 0 :=
      not_lt_of_ge hnonneg
    simp only [hremNonneg, ↓reduceIte]
    change ((coefficient.emod p.toNat).toNat.toUInt64).toNat < p.toNat
    rw [hword]
    exact hltNatP

theorem zpOfInt_toZMod (coefficient : Int) (p : UInt64)
    (hp : 0 < p.toNat) :
    CLPoly.Math.Zp.toZMod p.toNat (Zp.ofInt coefficient p) =
      (coefficient : ZMod p.toNat) := by
  unfold CLPoly.Math.Zp.toZMod Zp.ofInt
  have hpInt : (0 : Int) < p.toNat := by exact_mod_cast hp
  have hnonneg : 0 ≤ coefficient.emod p.toNat :=
    Int.emod_nonneg _ (by omega)
  have hltInt : coefficient.emod p.toNat < p.toNat :=
    Int.emod_lt_of_pos _ hpInt
  have hltNatP : (coefficient.emod p.toNat).toNat < p.toNat :=
    (Int.toNat_lt hnonneg).2 hltInt
  have hltSize : (coefficient.emod p.toNat).toNat < UInt64.size :=
    lt_trans hltNatP (UInt64.toNat_lt_size p)
  have hword : ((coefficient.emod p.toNat).toNat.toUInt64).toNat =
      (coefficient.emod p.toNat).toNat :=
    UInt64.toNat_ofNat_of_lt hltSize
  have hremNonneg : ¬coefficient % (p.toNat : Int) < 0 :=
    not_lt_of_ge hnonneg
  simp only [hremNonneg, ↓reduceIte]
  change (((coefficient.emod p.toNat).toNat.toUInt64).toNat :
      ZMod p.toNat) = (coefficient : ZMod p.toNat)
  rw [hword]
  have htoInt : ((coefficient.emod p.toNat).toNat : Int) =
      coefficient.emod p.toNat := Int.toNat_of_nonneg hnonneg
  calc
    ((coefficient.emod p.toNat).toNat : ZMod p.toNat) =
        (((coefficient.emod p.toNat).toNat : Int) : ZMod p.toNat) := by
          exact_mod_cast rfl
    _ = (coefficient.emod p.toNat : Int) := by rw [htoInt]
    _ = (coefficient : ZMod p.toNat) :=
      ZMod.intCast_mod coefficient p.toNat

/-- Pointwise effect of the source leading-coefficient baking block. -/
def HenselAdjustUnitRel (p : Nat)
    (before after : Array SparsePolyZp) : Prop :=
  before.size = after.size ∧
  ∀ index (hbefore : index < before.size) (hafter : index < after.size),
    ∃ scale : ZMod p, IsUnit scale ∧
      CLPoly.Math.SparsePolyZp.toPoly p after[index] =
        Polynomial.C scale *
          CLPoly.Math.SparsePolyZp.toPoly p before[index]

theorem HenselAdjustUnitRel.irreducible
    {p : Nat} {before after : Array SparsePolyZp}
    (hrel : HenselAdjustUnitRel p before after)
    (hirreducible : ∀ index (hindex : index < before.size),
      Irreducible (CLPoly.Math.SparsePolyZp.toPoly p before[index])) :
    ∀ index (hindex : index < after.size),
      Irreducible (CLPoly.Math.SparsePolyZp.toPoly p after[index]) := by
  intro index hafter
  have hbefore : index < before.size := by rw [hrel.1]; exact hafter
  rcases hrel.2 index hbefore hafter with ⟨scale, hscale, heq⟩
  rw [heq]
  exact (irreducible_isUnit_mul (Polynomial.isUnit_C.mpr hscale)).2
    (hirreducible index hbefore)

theorem HenselAdjustFirstFactorCorrect.unitRel
    {f : SparsePolyZZ} {factors adjusted : Array SparsePolyZp}
    {p : UInt64} [Fact (Nat.Prime p.toNat)]
    (hp2 : p.toNat * p.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical p.toNat factor)
    (hleadingNonzero : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod p.toNat) ≠ 0)
    (hcorrect : HenselAdjustFirstFactorCorrect f factors p adjusted) :
    HenselAdjustUnitRel p.toNat factors adjusted := by
  cases hcorrect with
  | adjusted leading first value hsource hfirst hvalue =>
      have hp : 0 < p.toNat := Nat.Prime.pos Fact.out
      have hzero : 0 < factors.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hfirst
        contradiction
      have hfirstGet : factors[0] = first :=
        Option.some.inj ((Array.getElem?_eq_getElem hzero).symm.trans hfirst)
      have hfirstMem : first ∈ factors.toList := by
        rw [← hfirstGet]
        exact Array.getElem_mem_toList hzero
      have hfirstCanonical := hfactors first hfirstMem
      let coefficient := Zp.ofInt leading.2 p
      have hcoefficient : CLPoly.Math.Zp.Reduced p.toNat coefficient :=
        zpOfInt_reduced leading.2 p hp
      have hcoefficientValue :
          CLPoly.Math.Zp.toZMod p.toNat coefficient =
            (leading.2 : ZMod p.toNat) :=
        zpOfInt_toZMod leading.2 p hp
      have hadjustedPoly :
          CLPoly.Math.SparsePolyZp.toPoly p.toNat value =
            Polynomial.C (leading.2 : ZMod p.toNat) *
              CLPoly.Math.SparsePolyZp.toPoly p.toNat first := by
        rw [hvalue, henselAdjustNormalization_toPoly,
          scaleZpCoeffs_toPoly p.toNat hp2 coefficient hcoefficient first
            hfirstCanonical.1, hcoefficientValue]
      have hscaleUnit : IsUnit (leading.2 : ZMod p.toNat) :=
        isUnit_iff_ne_zero.mpr (hleadingNonzero leading hsource)
      refine ⟨by simp [Array.set!, hzero], ?_⟩
      intro index hbefore hafter
      by_cases hindex : index = 0
      · subst index
        refine ⟨(leading.2 : ZMod p.toNat), hscaleUnit, ?_⟩
        simpa [Array.set!, hzero, hfirstGet] using hadjustedPoly
      · refine ⟨1, isUnit_one, ?_⟩
        have hzeroIndex : 0 ≠ index := by omega
        have hget : (factors.setIfInBounds 0 value)[index] =
            factors[index] := by
          simpa [hzeroIndex] using
            (@Array.getElem_setIfInBounds _ factors 0 value index hbefore)
        change CLPoly.Math.SparsePolyZp.toPoly p.toNat
            (factors.setIfInBounds 0 value)[index] = _
        rw [hget]
        simp

/-- Exact product effect of the source leading-coefficient baking block.
Unlike `unitRel`, this theorem retains the concrete scalar written by the
generated C++ execution and is therefore strong enough to seed the Hensel
tree product invariant without an associatedness oracle. -/
theorem HenselAdjustFirstFactorCorrect.product_eq
    {f : SparsePolyZZ} {factors adjusted : Array SparsePolyZp}
    {p : UInt64} [Fact (Nat.Prime p.toNat)]
    (hp2 : p.toNat * p.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical p.toNat factor)
    (hcorrect : HenselAdjustFirstFactorCorrect f factors p adjusted) :
    ∃ leading, f[0]? = some leading ∧
      (adjusted.toList.map
        (CLPoly.Math.SparsePolyZp.toPoly p.toNat)).prod =
        Polynomial.C (leading.2 : ZMod p.toNat) *
          (factors.toList.map
            (CLPoly.Math.SparsePolyZp.toPoly p.toNat)).prod := by
  cases hcorrect with
  | adjusted leading first value hsource hfirst hvalue =>
      have hp : 0 < p.toNat := Nat.Prime.pos Fact.out
      have hsize : 0 < factors.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hfirst
        contradiction
      have hfirstGet : factors[0] = first :=
        Option.some.inj ((Array.getElem?_eq_getElem hsize).symm.trans hfirst)
      have hfirstMem : first ∈ factors.toList := by
        rw [← hfirstGet]
        exact Array.getElem_mem_toList hsize
      have hfirstCanonical := hfactors first hfirstMem
      let coefficient := Zp.ofInt leading.2 p
      have hcoefficient : CLPoly.Math.Zp.Reduced p.toNat coefficient :=
        zpOfInt_reduced leading.2 p hp
      have hvaluePoly :
          CLPoly.Math.SparsePolyZp.toPoly p.toNat value =
            Polynomial.C (leading.2 : ZMod p.toNat) *
              CLPoly.Math.SparsePolyZp.toPoly p.toNat first := by
        rw [hvalue, henselAdjustNormalization_toPoly,
          scaleZpCoeffs_toPoly p.toNat hp2 coefficient hcoefficient first
            hfirstCanonical.1,
          zpOfInt_toZMod leading.2 p hp]
      have hlist : factors.toList = first :: factors.toList.drop 1 := by
        rw [← hfirstGet]
        simpa using (List.getElem_cons_drop (as := factors.toList) (i := 0)
          (by simpa using hsize)).symm
      refine ⟨leading, hsource, ?_⟩
      change (List.map (CLPoly.Math.SparsePolyZp.toPoly p.toNat)
        (factors.setIfInBounds 0 value).toList).prod = _
      rw [Array.toList_setIfInBounds, hlist]
      simp only [List.set_cons_zero, List.map_cons, List.prod_cons]
      rw [hvaluePoly]
      ring

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

/-- Semantic value selected by the two source target branches.  The default
branch names the concrete generated Mignotte computation, rather than an L2
oracle for the desired precision. -/
inductive HenselLiftTargetCorrect (f : SparsePolyZZ) (p : UInt64)
    (aTarget : Int32) : ZZ → Prop
  | mignotte (leading : UMonomial × ZZ)
      (hzero : aTarget = 0) (hleading : f[0]? = some leading) :
      HenselLiftTargetCorrect f p aTarget
        (2 * Int.natAbs leading.2 *
          (Generated.StrictHensel.__binomial_raw_ir
            (Int64.ofInt leading.1.deg)
            (Int64.ofInt leading.1.deg / 2) *
           Generated.StrictHensel.__isqrt_ceil_raw_ir
            (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f)))
  | explicit (hpositive : 0 < aTarget) :
      HenselLiftTargetCorrect f p aTarget
        ((p.toNat : ZZ) ^ aTarget.toNatClampNeg - 1)

/-- Raw-to-safe refinement for the complete target-selection conditional. -/
theorem __hensel_lift_target_raw_ir_refines
    (f : SparsePolyZZ) (p : UInt64) (aTarget : Int32)
    (hleading : ∃ leading, f[0]? = some leading)
    (hexponent : aTarget = 0 ∨ aTarget > 0) :
    ∃ target,
      Generated.StrictHensel.__hensel_lift_target_raw_ir f p aTarget =
        .ok target ∧
      HenselLiftTargetCorrect f p aTarget target := by
  rcases hleading with ⟨leading, hleading⟩
  by_cases hzero : aTarget = 0
  · subst aTarget
    rw [Generated.StrictHensel.__hensel_lift_target_raw_ir]
    simp only [ite_true, Generated.StrictHensel.__mignotte_bound_upoly_raw_ir,
      hleading, bind, Except.bind]
    exact ⟨_, rfl, .mignotte leading rfl hleading⟩
  · have hpositive : aTarget > 0 := hexponent.resolve_left hzero
    rcases __hensel_explicit_target_raw_ir_refines p aTarget hpositive with
      ⟨target, hrun, hcorrect⟩
    refine ⟨target, ?_, ?_⟩
    · simp [Generated.StrictHensel.__hensel_lift_target_raw_ir, hzero, hrun]
    · rw [HenselExplicitTargetCorrect] at hcorrect
      subst target
      exact .explicit hpositive

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
  inverseReady : ∀ state leading,
    HenselEEAPrefix ops initial state → state.r1.isEmpty = true →
    state.r0[0]? = some leading →
    ∃ inverse, ops.inverse leading.2 = .ok inverse
  scaleReady : ∀ state leading inverse,
    HenselEEAPrefix ops initial state → state.r1.isEmpty = true →
    state.r0[0]? = some leading →
    ops.inverse leading.2 = .ok inverse →
    ∃ gcd s t,
      ops.scaleNormalize inverse state.r0 = .ok gcd ∧
      ops.scaleNormalize inverse state.s0 = .ok s ∧
      ops.scaleNormalize inverse state.t0 = .ok t

/-- Exact semantic execution trace of the source extended-Euclidean loop. -/
inductive HenselEEACorrect
    (ops : Generated.StrictHensel.HenselEEARawOps) :
    Generated.StrictHensel.HenselEEAState →
      (SparsePolyZp × SparsePolyZp × SparsePolyZp) → Prop
  | done (state : Generated.StrictHensel.HenselEEAState)
      (leading : UMonomial × Zp) (inverse : Zp)
      (gcd s t : SparsePolyZp)
      (hdone : state.r1.isEmpty = true)
      (hleading : state.r0[0]? = some leading)
      (hinverse : ops.inverse leading.2 = .ok inverse)
      (hgcd : ops.scaleNormalize inverse state.r0 = .ok gcd)
      (hs : ops.scaleNormalize inverse state.s0 = .ok s)
      (ht : ops.scaleNormalize inverse state.t0 = .ok t) :
      HenselEEACorrect ops state (gcd, s, t)
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
    rcases hinvariant.inverseReady state leading hprefix hdone hleading with
      ⟨inverse, hinverse⟩
    rcases hinvariant.scaleReady state leading inverse hprefix hdone hleading
        hinverse with ⟨gcd, s, t, hgcd, hs, ht⟩
    rw [hleading]
    simp only
    rw [hinverse]
    simp only [hgcd, hs, ht, bind, Except.bind]
    exact ⟨(gcd, s, t), rfl,
      .done state leading inverse gcd s t hdone hleading hinverse hgcd hs ht⟩
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

def strictHenselEEAMulCoefficientIR (this : DenseUPolyZp)
    (left right : Zp) : Zp :=
  { val := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
      left.val right.val,
    prime := this._p }

def strictHenselEEAScaleNormalizeIR (this : DenseUPolyZp)
    (coefficient : Zp) (f : SparsePolyZp) : RawExec SparsePolyZp :=
  .ok (SparsePolyZp.normalization (f.map fun term =>
    (term.1, strictHenselEEAMulCoefficientIR this term.2 coefficient)))

def strictHenselEEARawOps (this : DenseUPolyZp) :
    Generated.StrictHensel.HenselEEARawOps where
  divmod := henselDivmodVHCIR this
  inverse := fun coefficient => .ok
      { val := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
          coefficient.val,
        prime := this._p }
  scaleNormalize := strictHenselEEAScaleNormalizeIR this

/-- The terminal EEA inverse call is the generated C++ `inv_prime` path, and
its returned runtime coefficient is the field inverse of the concrete
nonzero leading coefficient. -/
theorem strictHenselEEARawOps_inverse_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (coefficient : Zp)
    (hreduced : CLPoly.Math.Zp.Reduced this._p.toNat coefficient)
    (hnonzero : coefficient.val ≠ 0) :
    ∃ inverse,
      (strictHenselEEARawOps this).inverse coefficient = .ok inverse ∧
      CLPoly.Math.Zp.Reduced this._p.toNat inverse ∧
      CLPoly.Math.Zp.toZMod this._p.toNat inverse *
          CLPoly.Math.Zp.toZMod this._p.toNat coefficient = 1 := by
  let inverse : Zp :=
    { val := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        coefficient.val,
      prime := this._p }
  have hpositive : 0 < coefficient.val.toNat := by
    exact Nat.pos_of_ne_zero (fun hzero => hnonzero
      (UInt64.toNat_inj.mp (by simpa using hzero)))
  have hcorrect :=
    CLPoly.Impl.StrictWordArithmetic.dense_upoly_zp_nmod_inv_ir_correct this
      coefficient.val (Fact.out : Nat.Prime this._p.toNat) hpositive
      hreduced.2
  dsimp only at hcorrect
  refine ⟨inverse, rfl, ⟨rfl, hcorrect.1⟩, ?_⟩
  simpa [inverse, CLPoly.Math.Zp.toZMod, hreduced.1] using hcorrect.2

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
  inverseReady := by
    intro state leading _ _ _
    let inverse : Zp :=
      { val := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
          leading.2.val,
        prime := this._p }
    exact ⟨inverse, rfl⟩
  scaleReady := by
    intro state leading inverse _ _ _ _
    exact ⟨SparsePolyZp.normalization (state.r0.map fun term =>
        (term.1,
          { val := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
              term.2.val inverse.val,
            prime := this._p })),
      SparsePolyZp.normalization (state.s0.map fun term =>
        (term.1,
          { val := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
              term.2.val inverse.val,
            prime := this._p })),
      SparsePolyZp.normalization (state.t0.map fun term =>
        (term.1,
          { val := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
              term.2.val inverse.val,
            prime := this._p })), rfl, rfl, rfl⟩

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

theorem strictHenselEEAMulCoefficientIR_refines
    (this : DenseUPolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (left right : Zp)
    (hleft : CLPoly.Math.Zp.Reduced this._p.toNat left) :
    CLPoly.Math.Zp.Reduced this._p.toNat
        (strictHenselEEAMulCoefficientIR this left right) ∧
      CLPoly.Math.Zp.toZMod this._p.toNat
          (strictHenselEEAMulCoefficientIR this left right) =
        CLPoly.Math.Zp.toZMod this._p.toNat left *
          CLPoly.Math.Zp.toZMod this._p.toNat right := by
  have hp : 0 < this._p.toNat := by
    exact lt_of_le_of_lt (Nat.zero_le _) hleft.2
  have hmul :=
    CLPoly.Impl.StrictWordArithmetic.nmod_mul_ir_correct_of_configured this
      left.val right.val hcfg hleft.2
  have hbound :
      (Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this left.val
          right.val).toNat < this._p.toNat := by
    rw [hmul]
    exact Nat.mod_lt _ hp
  refine ⟨⟨rfl, hbound⟩, ?_⟩
  unfold strictHenselEEAMulCoefficientIR CLPoly.Math.Zp.toZMod
  dsimp only
  rw [hmul, ZMod.natCast_mod]
  push_cast
  rfl

private theorem strictHenselEEAScaleList_toPoly
    (this : DenseUPolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (coefficient : Zp) :
    ∀ terms : List (UMonomial × Zp),
      CLPoly.Math.SparsePolyZp.AllReduced this._p.toNat terms →
      CLPoly.Math.listSum this._p.toNat
          (terms.map fun term =>
            (term.1, strictHenselEEAMulCoefficientIR this term.2 coefficient)) =
        Polynomial.C (CLPoly.Math.Zp.toZMod this._p.toNat coefficient) *
          CLPoly.Math.listSum this._p.toNat terms := by
  intro terms hterms
  induction terms with
  | nil => simp [CLPoly.Math.listSum]
  | cons term rest ih =>
      have hterm := hterms term List.mem_cons_self
      have hrest : CLPoly.Math.SparsePolyZp.AllReduced this._p.toNat rest :=
        fun item hitem => hterms item (List.mem_cons_of_mem term hitem)
      have hcoefficient :=
        strictHenselEEAMulCoefficientIR_refines this hcfg term.2 coefficient
          hterm
      rcases term with ⟨monomial, value⟩
      simp only [List.map_cons, CLPoly.Math.listSum_cons]
      rw [hcoefficient.2, ih hrest]
      rw [mul_add, Polynomial.C_mul_monomial]
      ring

/-- Exact semantic contract for one of the three terminal source scaling
loops followed by sparse normalization. -/
theorem strictHenselEEAScaleNormalizeIR_refines
    (this : DenseUPolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (coefficient : Zp) (f : SparsePolyZp)
    (hf : CLPoly.Math.SparsePolyZp.WellFormed_arr this._p.toNat f) :
    ∃ output,
      strictHenselEEAScaleNormalizeIR this coefficient f = .ok output ∧
      CLPoly.Math.SparsePolyZp.WellFormed_arr this._p.toNat output ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output =
        Polynomial.C (CLPoly.Math.Zp.toZMod this._p.toNat coefficient) *
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat f := by
  let mapped := f.map fun term =>
    (term.1, strictHenselEEAMulCoefficientIR this term.2 coefficient)
  let output := SparsePolyZp.normalization mapped
  have hmapped : CLPoly.Math.SparsePolyZp.WellFormed_arr this._p.toNat
      mapped := by
    intro term hterm
    dsimp only [mapped] at hterm
    rw [Array.toList_map, List.mem_map] at hterm
    rcases hterm with ⟨source, hsource, rfl⟩
    exact (strictHenselEEAMulCoefficientIR_refines this hcfg source.2
      coefficient (hf source hsource)).1
  refine ⟨output, rfl,
    henselEEANormalization_wellFormed this._p.toNat mapped hmapped, ?_⟩
  rw [henselEEANormalization_toPoly this._p.toNat mapped]
  unfold CLPoly.Math.SparsePolyZp.toPoly
  rw [show mapped.toList = f.toList.map (fun term =>
      (term.1, strictHenselEEAMulCoefficientIR this term.2 coefficient)) by
    simp [mapped]]
  exact strictHenselEEAScaleList_toPoly this hcfg coefficient f.toList hf

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

/-- The three arrays returned by the concrete generated EEA satisfy the
Bézout equation.  The proof follows the actual successful raw trace; at the
terminal state it uses the three generated coefficient-scaling executions. -/
theorem strictHenselEEACorrect_bezout
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (left right : SparsePolyZp)
    (state : Generated.StrictHensel.HenselEEAState)
    (output : SparsePolyZp × SparsePolyZp × SparsePolyZp)
    (hstate : StrictHenselEEAStateInvariant this._p.toNat state)
    (halgebra : StrictHenselEEAAlgebraicInvariant this._p.toNat left right
      state)
    (hcorrect : HenselEEACorrect (strictHenselEEARawOps this) state output) :
    CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 =
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.2.1 *
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left +
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.2.2 *
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat right := by
  induction hcorrect with
  | done state leading inverse gcd s t hdone hleading hinverse hgcd hs ht =>
      rcases strictHenselEEAScaleNormalizeIR_refines this hcfg inverse
          state.r0 hstate.r0Canonical.1 with
        ⟨gcd', hgcdRun, _, hgcdPoly⟩
      rcases strictHenselEEAScaleNormalizeIR_refines this hcfg inverse
          state.s0 halgebra.s0WellFormed with
        ⟨s', hsRun, _, hsPoly⟩
      rcases strictHenselEEAScaleNormalizeIR_refines this hcfg inverse
          state.t0 halgebra.t0WellFormed with
        ⟨t', htRun, _, htPoly⟩
      have hgcdEq : gcd' = gcd := Except.ok.inj (hgcdRun.symm.trans (by
        simpa [strictHenselEEARawOps] using hgcd))
      have hsEq : s' = s := Except.ok.inj (hsRun.symm.trans (by
        simpa [strictHenselEEARawOps] using hs))
      have htEq : t' = t := Except.ok.inj (htRun.symm.trans (by
        simpa [strictHenselEEARawOps] using ht))
      subst gcd'
      subst s'
      subst t'
      rw [hgcdPoly, hsPoly, htPoly, halgebra.r0Equation]
      ring
  | step state quotient remainder output hcontinue hrun htail ih =>
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
        ⟨actual, hactualRun, hactualCorrect⟩
      have hrun' : henselDivmodVHCIR this state.r0 state.r1 =
          .ok (quotient, remainder) := by
        simpa [strictHenselEEARawOps] using hrun
      have hactual : actual = (quotient, remainder) :=
        Except.ok.inj (hactualRun.symm.trans hrun')
      subst actual
      have hnextState : StrictHenselEEAStateInvariant this._p.toNat
          (Generated.StrictHensel.henselEEANextState state quotient
            remainder) :=
        ⟨hstate.r1Canonical, hactualCorrect.remainderCanonical, hr1⟩
      have hnextAlgebra := StrictHenselEEAAlgebraicInvariant.step
        this._p.toNat left right state quotient remainder h2p hp2 hr1
        halgebra hactualCorrect
      exact ih hnextState hnextAlgebra

theorem henselEEAToPoly_leadingCoeff_eq_head
    (p : Nat) [Fact (Nat.Prime p)] (f : SparsePolyZp)
    (hcanonical : CLPoly.Math.SparsePolyZp.Canonical p f)
    (hnonempty : 0 < f.size) :
    (CLPoly.Math.SparsePolyZp.toPoly p f).leadingCoeff =
      CLPoly.Math.Zp.toZMod p f[0].2 := by
  have hdegree :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      p f hcanonical hnonempty
  have hnatDegree : (CLPoly.Math.SparsePolyZp.toPoly p f).natDegree =
      f[0].1.deg := Polynomial.natDegree_eq_of_degree_eq_some hdegree
  have hlistNonempty : f.toList ≠ [] := by
    intro hempty
    have hlength := congrArg List.length hempty
    have hsizeZero : f.size = 0 := by simpa using hlength
    exact (Nat.ne_of_gt hnonempty) hsizeZero
  obtain ⟨head, rest, hlist⟩ := List.exists_cons_of_ne_nil hlistNonempty
  have hheadEq : head = f[0] := by
    have hget := Array.getElem_toList hnonempty
    simpa [hlist] using hget
  have hchain : List.IsChain
      (fun a b : UMonomial × Zp => a.1.deg > b.1.deg)
      (head :: rest) := by
    simpa [hlist] using hcanonical.2.1
  rw [Polynomial.leadingCoeff, hnatDegree]
  unfold CLPoly.Math.SparsePolyZp.toPoly
  rw [hlist, ← hheadEq]
  exact CLPoly.Math.listSum_coeff_at_head p head rest hchain

/-- The first component returned by the generated EEA is monic.  Its terminal
leading coefficient is normalized by the actual generated inverse and
multiplication calls. -/
theorem strictHenselEEACorrect_monic
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (state : Generated.StrictHensel.HenselEEAState)
    (output : SparsePolyZp × SparsePolyZp × SparsePolyZp)
    (hstate : StrictHenselEEAStateInvariant this._p.toNat state)
    (hcorrect : HenselEEACorrect (strictHenselEEARawOps this) state output) :
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1).Monic := by
  induction hcorrect with
  | done state leading inverse gcd s t hdone hleading hinverse hgcd hs ht =>
      let head := state.r0[0]'hstate.r0Nonempty
      have hheadOpt : state.r0[0]? = some head := by
        simpa [head] using Array.getElem?_eq_getElem hstate.r0Nonempty
      have hleadingEq : leading = head :=
        Option.some.inj (hleading.symm.trans hheadOpt)
      subst leading
      have hheadMem : head ∈ state.r0.toList := by
        simpa [head] using Array.getElem_mem_toList state.r0 0
          hstate.r0Nonempty
      have hheadReduced := hstate.r0Canonical.1 head hheadMem
      have hheadNonzero := hstate.r0Canonical.2.2 head hheadMem
      rcases strictHenselEEARawOps_inverse_refines this head.2
          hheadReduced hheadNonzero with
        ⟨inverse', hinverseRun, _, hinverseField⟩
      have hinverseEq : inverse' = inverse := Except.ok.inj
        (hinverseRun.symm.trans hinverse)
      subst inverse'
      rcases strictHenselEEAScaleNormalizeIR_refines this hcfg inverse
          state.r0 hstate.r0Canonical.1 with
        ⟨gcd', hgcdRun, _, hgcdPoly⟩
      have hgcdEq : gcd' = gcd := Except.ok.inj (hgcdRun.symm.trans (by
        simpa [strictHenselEEARawOps] using hgcd))
      subst gcd'
      rw [hgcdPoly]
      apply Polynomial.monic_C_mul_of_mul_leadingCoeff_eq_one
      rw [henselEEAToPoly_leadingCoeff_eq_head this._p.toNat state.r0
        hstate.r0Canonical hstate.r0Nonempty]
      exact hinverseField
  | step state quotient remainder output hcontinue hrun htail ih =>
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
        ⟨actual, hactualRun, hactualCorrect⟩
      have hrun' : henselDivmodVHCIR this state.r0 state.r1 =
          .ok (quotient, remainder) := by
        simpa [strictHenselEEARawOps] using hrun
      have hactual : actual = (quotient, remainder) :=
        Except.ok.inj (hactualRun.symm.trans hrun')
      subst actual
      exact ih ⟨hstate.r1Canonical, hactualCorrect.remainderCanonical, hr1⟩

/-- The concrete first EEA result divides both remainder registers of the
state from which its raw trace started.  The step case reconstructs the old
dividend from the actual generated quotient/remainder equation. -/
theorem strictHenselEEACorrect_dvd
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (state : Generated.StrictHensel.HenselEEAState)
    (output : SparsePolyZp × SparsePolyZp × SparsePolyZp)
    (hstate : StrictHenselEEAStateInvariant this._p.toNat state)
    (hcorrect : HenselEEACorrect (strictHenselEEARawOps this) state output) :
    CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 ∣
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r0 ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 ∣
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r1 := by
  induction hcorrect with
  | done state leading inverse gcd s t hdone hleading hinverse hgcd hs ht =>
      let head := state.r0[0]'hstate.r0Nonempty
      have hheadOpt : state.r0[0]? = some head := by
        simpa [head] using Array.getElem?_eq_getElem hstate.r0Nonempty
      have hleadingEq : leading = head :=
        Option.some.inj (hleading.symm.trans hheadOpt)
      subst leading
      have hheadMem : head ∈ state.r0.toList := by
        simpa [head] using Array.getElem_mem_toList state.r0 0
          hstate.r0Nonempty
      have hheadReduced := hstate.r0Canonical.1 head hheadMem
      have hheadNonzero := hstate.r0Canonical.2.2 head hheadMem
      rcases strictHenselEEARawOps_inverse_refines this head.2
          hheadReduced hheadNonzero with
        ⟨inverse', hinverseRun, _, hinverseField⟩
      have hinverseEq : inverse' = inverse := Except.ok.inj
        (hinverseRun.symm.trans hinverse)
      subst inverse'
      rcases strictHenselEEAScaleNormalizeIR_refines this hcfg inverse
          state.r0 hstate.r0Canonical.1 with
        ⟨gcd', hgcdRun, _, hgcdPoly⟩
      have hgcdEq : gcd' = gcd := Except.ok.inj (hgcdRun.symm.trans (by
        simpa [strictHenselEEARawOps] using hgcd))
      subst gcd'
      constructor
      · refine ⟨Polynomial.C (CLPoly.Math.Zp.toZMod this._p.toNat head.2),
          ?_⟩
        rw [hgcdPoly]
        calc
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r0 =
              1 * CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r0 := by
                ring
          _ = Polynomial.C
                (CLPoly.Math.Zp.toZMod this._p.toNat inverse *
                  CLPoly.Math.Zp.toZMod this._p.toNat head.2) *
                CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r0 := by
              rw [hinverseField, Polynomial.C_1]
          _ = (Polynomial.C (CLPoly.Math.Zp.toZMod this._p.toNat inverse) *
                CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r0) *
              Polynomial.C
                (CLPoly.Math.Zp.toZMod this._p.toNat head.2) := by
              rw [Polynomial.C_mul]
              ring
      · have hr1Empty : state.r1 = #[] := by
          simpa [Array.isEmpty_iff] using hdone
        rw [hr1Empty, CLPoly.Math.SparsePolyZp.toPoly_empty]
        exact dvd_zero _
  | step state quotient remainder output hcontinue hrun htail ih =>
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
        ⟨actual, hactualRun, hactualCorrect⟩
      have hrun' : henselDivmodVHCIR this state.r0 state.r1 =
          .ok (quotient, remainder) := by
        simpa [strictHenselEEARawOps] using hrun
      have hactual : actual = (quotient, remainder) :=
        Except.ok.inj (hactualRun.symm.trans hrun')
      subst actual
      have htailDvd := ih
        ⟨hstate.r1Canonical, hactualCorrect.remainderCanonical, hr1⟩
      rcases htailDvd.1 with ⟨rightWitness, hrightWitness⟩
      rcases htailDvd.2 with ⟨remainderWitness, hremainderWitness⟩
      change CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.r1 =
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 *
          rightWitness at hrightWitness
      change CLPoly.Math.SparsePolyZp.toPoly this._p.toNat remainder =
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 *
          remainderWitness at hremainderWitness
      constructor
      · refine ⟨CLPoly.Math.SparsePolyZp.toPoly this._p.toNat quotient *
            rightWitness + remainderWitness, ?_⟩
        rw [← hactualCorrect.equation, hrightWitness, hremainderWitness]
        ring
      · exact ⟨rightWitness, hrightWitness⟩

theorem henselEEAOne_refines (p : UInt64) [Fact (Nat.Prime p.toNat)] :
    CLPoly.Math.SparsePolyZp.Canonical p.toNat
        (#[(UMonomial.mk 0, Zp.ofUInt64 1 p)] : SparsePolyZp) ∧
      CLPoly.Math.SparsePolyZp.toPoly p.toNat
        (#[(UMonomial.mk 0, Zp.ofUInt64 1 p)] : SparsePolyZp) = 1 := by
  have hp : 1 < p.toNat := (Fact.out : Nat.Prime p.toNat).one_lt
  have hpWord : p ≠ 0 := by
    intro hzero
    subst p
    norm_num at hp
  have hmodNat : 1 % p.toNat = 1 := Nat.mod_eq_of_lt hp
  have hpWordLt : (1 : UInt64) < p := by exact_mod_cast hp
  have hmodWord : (1 : UInt64) % p = 1 := by
    exact UInt64.mod_eq_of_lt hpWordLt
  simp [CLPoly.Math.SparsePolyZp.Canonical,
    CLPoly.Math.SparsePolyZp.WellFormed_arr,
    CLPoly.Math.SparsePolyZp.AllReduced,
    CLPoly.Math.SparsePolyZp.toPoly, CLPoly.Math.listSum,
    CLPoly.Math.Zp.Reduced, Zp.ofUInt64, CLPoly.Math.Zp.toZMod,
    hp, hmodNat, hmodWord]

/-- The six source initial assignments establish both execution safety and
the two base Bézout identities. -/
theorem strictHenselEEAInitialState_invariants
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (left right : SparsePolyZp)
    (hleftCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat right)
    (hleftNonempty : 0 < left.size) :
    StrictHenselEEAStateInvariant this._p.toNat
        (Generated.StrictHensel.henselEEAInitialState this._p left right) ∧
      StrictHenselEEAAlgebraicInvariant this._p.toNat left right
        (Generated.StrictHensel.henselEEAInitialState this._p left right) := by
  have hone := henselEEAOne_refines this._p
  constructor
  · exact ⟨hleftCanonical, hrightCanonical, hleftNonempty⟩
  · dsimp [Generated.StrictHensel.henselEEAInitialState]
    refine ⟨hone.1.1, ?_, ?_, hone.1.1, ?_, ?_⟩
    · intro term hterm
      dsimp [Generated.StrictHensel.henselEEAInitialState] at hterm
      simp at hterm
    · intro term hterm
      dsimp [Generated.StrictHensel.henselEEAInitialState] at hterm
      simp at hterm
    · simp only [Generated.StrictHensel.HenselEEAState.r0,
        Generated.StrictHensel.HenselEEAState.s0,
        Generated.StrictHensel.HenselEEAState.t0,
        CLPoly.Math.SparsePolyZp.toPoly_empty, zero_mul, add_zero]
      rw [hone.2]
      simp
    · simp only [Generated.StrictHensel.HenselEEAState.r1,
        Generated.StrictHensel.HenselEEAState.s1,
        Generated.StrictHensel.HenselEEAState.t1,
        CLPoly.Math.SparsePolyZp.toPoly_empty, zero_mul, zero_add]
      rw [hone.2]
      simp

theorem polynomial_eq_gcd_of_monic_bezout_dvd
    {p : Nat} [Fact (Nat.Prime p)]
    (left right gcd s t : Polynomial (ZMod p))
    (hleftNonzero : left ≠ 0)
    (hbezout : gcd = s * left + t * right)
    (hmonic : gcd.Monic)
    (hdvdLeft : gcd ∣ left) (hdvdRight : gcd ∣ right) :
    gcd = GCDMonoid.gcd left right := by
  let standard := GCDMonoid.gcd left right
  have hgcdDvd : gcd ∣ standard := by
    exact GCDMonoid.dvd_gcd hdvdLeft hdvdRight
  have hstandardDvd : standard ∣ gcd := by
    rcases GCDMonoid.gcd_dvd_left left right with ⟨leftWitness, hleft⟩
    rcases GCDMonoid.gcd_dvd_right left right with ⟨rightWitness, hright⟩
    change left = standard * leftWitness at hleft
    change right = standard * rightWitness at hright
    refine ⟨s * leftWitness + t * rightWitness, ?_⟩
    calc
      gcd = s * left + t * right := hbezout
      _ = s * (standard * leftWitness) +
          t * (standard * rightWitness) :=
        congrArg₂ (fun a b => s * a + t * b) hleft hright
      _ = standard * (s * leftWitness + t * rightWitness) := by ring
  have hstandardNonzero : standard ≠ 0 := by
    intro hzero
    have hdvd := GCDMonoid.gcd_dvd_left left right
    change standard ∣ left at hdvd
    rw [hzero, zero_dvd_iff] at hdvd
    exact hleftNonzero hdvd
  have hstandardMonic : standard.Monic := by
    change (GCDMonoid.gcd left right).Monic
    rw [← normalize_gcd left right]
    exact Polynomial.monic_normalize hstandardNonzero
  change gcd = standard
  exact Polynomial.eq_of_monic_of_associated hmonic hstandardMonic
    (associated_of_dvd_dvd hgcdDvd hstandardDvd)

/-- Entry-level semantic contract for the generated, well-founded EEA.  The
returned concrete arrays are produced by the raw program and satisfy its
Bézout identity over `ZMod p`. -/
theorem strictHenselEEAEntryIR_refines_gcd
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (left right : SparsePolyZp)
    (hleftCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat right)
    (hleftNonempty : 0 < left.size) :
    ∃ gcd s t,
      Generated.StrictHensel.__polynomial_GCD_eea_raw_ir
          (strictHenselEEARawOps this) (strictHenselEEATermination this)
          (Generated.StrictHensel.henselEEAInitialState this._p left right) =
        .ok (gcd, s, t) ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd =
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat s *
            CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left +
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat t *
            CLPoly.Math.SparsePolyZp.toPoly this._p.toNat right ∧
      (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd).Monic ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd ∣
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd ∣
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat right ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd =
        GCDMonoid.gcd
          (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left)
          (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat right) := by
  let initial :=
    Generated.StrictHensel.henselEEAInitialState this._p left right
  have hinitial := strictHenselEEAInitialState_invariants this left right
    hleftCanonical hrightCanonical hleftNonempty
  rcases strictHenselEEAIR_refines this hcfg initial hinitial.1 with
    ⟨⟨gcd, s, t⟩, hrun, hcorrect⟩
  have hdvd := strictHenselEEACorrect_dvd this hcfg initial (gcd, s, t)
    hinitial.1 hcorrect
  have hbezout := strictHenselEEACorrect_bezout this hcfg h2p hp2 left right initial
      (gcd, s, t) hinitial.1 hinitial.2 hcorrect
  have hmonic := strictHenselEEACorrect_monic this hcfg initial (gcd, s, t)
      hinitial.1 hcorrect
  have hdvdLeft : CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd ∣
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left := by
    simpa [initial, Generated.StrictHensel.henselEEAInitialState] using hdvd.1
  have hdvdRight : CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd ∣
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat right := by
    simpa [initial, Generated.StrictHensel.henselEEAInitialState] using hdvd.2
  have hleftNonzero :
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left ≠ 0 := by
    intro hzero
    have hdegree :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
        this._p.toNat left hleftCanonical hleftNonempty
    rw [hzero, Polynomial.degree_zero] at hdegree
    simp at hdegree
  refine ⟨gcd, s, t, hrun, hbezout, hmonic, hdvdLeft, hdvdRight, ?_⟩
  exact polynomial_eq_gcd_of_monic_bezout_dvd
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat left)
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat right)
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat gcd)
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat s)
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat t)
    hleftNonzero hbezout hmonic hdvdLeft hdvdRight

def strictHenselTreeBuildRawOps (this : DenseUPolyZp)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this) :
    Generated.StrictHensel.HenselTreeBuildRawOps where
  mul := fun left right => StrictDDF.strictMulIR this left right mulProvider
  eea := fun left right =>
    Generated.StrictHensel.__polynomial_GCD_eea_raw_ir
      (strictHenselEEARawOps this) (strictHenselEEATermination this)
      (Generated.StrictHensel.henselEEAInitialState this._p left right)

/-- The source `poly_convert` representation change preserves the polynomial
over the input prime.  This is the bridge from the concrete integer-valued
tree fields back to the ZMod semantics used by the EEA proof. -/
theorem henselTreeZpToZZIR_toPolyMod (p : Nat) (f : SparsePolyZp) :
    toPolyMod p (Generated.StrictHensel.henselTreeZpToZZIR f) =
      CLPoly.Math.SparsePolyZp.toPoly p f := by
  unfold toPolyMod Generated.StrictHensel.henselTreeZpToZZIR
    SparsePolyZZ.toPoly CLPoly.Math.SparsePolyZp.toPoly
  rw [Array.toList_map]
  induction f.toList with
  | nil => simp [CLPoly.Math.listSum]
  | cons term rest ih =>
      rcases term with ⟨monomial, coefficient⟩
      rw [List.map_map]
      simp only [List.map_cons, List.sum_cons, map_add,
        Polynomial.map_monomial, CLPoly.Math.listSum_cons]
      have ih' : Polynomial.map (Int.castRingHom (ZMod p))
          (rest.map fun term =>
            Polynomial.monomial term.1.deg (term.2.val.toNat : Int)).sum =
          CLPoly.Math.listSum p rest := by
        simpa only [List.map_map, Function.comp_apply] using ih
      have htail := congrArg
        (fun tail => Polynomial.monomial monomial.deg
          (CLPoly.Math.Zp.toZMod p coefficient) + tail) ih'
      simpa [List.map_map, Function.comp_def, CLPoly.Math.Zp.toZMod] using htail

/-- The generated coefficient representation conversion used by the C++ tree
builder preserves strict degree order and nonzero coefficients. -/
theorem henselTreeZpToZZIR_canonical (p : Nat) (f : SparsePolyZp)
    (hf : CLPoly.Math.SparsePolyZp.Canonical p f) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (Generated.StrictHensel.henselTreeZpToZZIR f) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  constructor
  · rw [Generated.StrictHensel.henselTreeZpToZZIR, Array.toList_map,
      List.isChain_map]
    simpa using hf.2.1
  · intro output houtput
    rw [Generated.StrictHensel.henselTreeZpToZZIR, Array.toList_map,
      List.mem_map] at houtput
    rcases houtput with ⟨source, hsource, rfl⟩
    have hnonzero := hf.2.2 source hsource
    simp only [Prod.snd]
    apply Int.ofNat_ne_zero.mpr
    intro hzero
    exact hnonzero (UInt64.toNat_inj.mp (by simp [hzero]))

/-- A canonical nonempty finite-field polynomial whose decoded polynomial is
monic converts to an integer sparse array with a literal coefficient-one
head.  This is a representation theorem about the generated `poly_convert`,
not merely a modular polynomial equality. -/
theorem henselTreeZpToZZIR_hasPhysicalOneHead_of_monic
    (p : Nat) [Fact (Nat.Prime p)] (f : SparsePolyZp)
    (hf : CLPoly.Math.SparsePolyZp.Canonical p f)
    (hnonempty : 0 < f.size)
    (hmonic : (CLPoly.Math.SparsePolyZp.toPoly p f).Monic) :
    HasPhysicalOneHead
      (Generated.StrictHensel.henselTreeZpToZZIR f) := by
  have hleading := henselEEAToPoly_leadingCoeff_eq_head p f hf hnonempty
  rw [hmonic.leadingCoeff] at hleading
  have hheadMem : f[0] ∈ f.toList := Array.getElem_mem_toList hnonempty
  have hheadReduced := hf.1 f[0] hheadMem
  have hvalNat : f[0].2.val.toNat = 1 := by
    have hval := congrArg ZMod.val hleading
    unfold CLPoly.Math.Zp.toZMod at hval
    rw [ZMod.val_natCast, Nat.mod_eq_of_lt hheadReduced.2] at hval
    have hp : 1 < p := (Fact.out : Nat.Prime p).one_lt
    simpa [ZMod.val_one, Nat.mod_eq_of_lt hp] using hval.symm
  have hlist : f.toList = f[0] :: f.toList.drop 1 := by
    simpa using (List.getElem_cons_drop (as := f.toList) (i := 0)
      (by simpa using hnonempty)).symm
  refine ⟨f[0].1,
    (f.toList.drop 1).map fun term =>
      (term.1, (term.2.val.toNat : Int)), ?_⟩
  rw [Generated.StrictHensel.henselTreeZpToZZIR, Array.toList_map, hlist,
    List.map_cons]
  simp [hvalNat]

/-- Mathematical denotation of the same half-open source interval consumed
by a tree-product loop.  This is used only in the proof layer. -/
noncomputable def henselFactorRangeProduct (p : Nat)
    (factors : Array SparsePolyZp) (stop : Nat) : Nat →
      Polynomial (ZMod p)
  | index =>
      if index < stop then
        CLPoly.Math.SparsePolyZp.toPoly p factors[index]! *
          henselFactorRangeProduct p factors stop (index + 1)
      else 1
termination_by index => stop - index
decreasing_by simp_wf; omega

/-- A half-open product of the actual factor array is monic when every
factor read by that interval is monic. -/
theorem henselFactorRangeProduct_monic
    (p : Nat) (factors : Array SparsePolyZp) (start stop : Nat)
    (hstop : stop ≤ factors.size)
    (hmonic : ∀ index (hindex : index < factors.size),
      start ≤ index → index < stop →
      (CLPoly.Math.SparsePolyZp.toPoly p
        (getElem factors index hindex)).Monic) :
    (henselFactorRangeProduct p factors stop start).Monic := by
  rw [henselFactorRangeProduct]
  by_cases hmore : start < stop
  · rw [if_pos hmore]
    have hstart : start < factors.size := lt_of_lt_of_le hmore hstop
    rw [getElem!_pos factors start hstart]
    exact (hmonic start hstart (by omega) hmore).mul
      (henselFactorRangeProduct_monic p factors (start + 1) stop hstop
        (by
          intro index hindex hstartIndex hindexStop
          exact hmonic index hindex (by omega) hindexStop))
  · rw [if_neg hmore]
    exact Polynomial.monic_one
termination_by stop - start
decreasing_by simp_wf; omega

/-- Ordered list denoted by the same half-open factor interval as
`henselFactorRangeProduct`.  This retains the source array order and is used
to state the pointwise leaf-origin theorem; it does not predict any lifted
factor values. -/
def henselFactorRangeList (factors : Array SparsePolyZp) (stop : Nat) :
    Nat → List SparsePolyZp
  | index =>
      if index < stop then
        factors[index]! :: henselFactorRangeList factors stop (index + 1)
      else []
termination_by index => stop - index
decreasing_by simp_wf; omega

@[simp] theorem henselFactorRangeList_self
    (factors : Array SparsePolyZp) (stop : Nat) :
    henselFactorRangeList factors stop stop = [] := by
  rw [henselFactorRangeList]
  simp

/-- The source interval ending at the array size is its concrete list suffix. -/
theorem henselFactorRangeList_suffix (factors : Array SparsePolyZp)
    (index : Nat) :
    henselFactorRangeList factors factors.size index =
      factors.toList.drop index := by
  rw [henselFactorRangeList]
  by_cases hindex : index < factors.size
  · rw [if_pos hindex]
    have hdrop := List.drop_eq_getElem_cons
      (l := factors.toList) (i := index) (by simpa using hindex)
    rw [henselFactorRangeList_suffix factors (index + 1)]
    simpa [getElem!_pos factors index hindex] using hdrop.symm
  · rw [if_neg hindex]
    exact (List.drop_eq_nil_iff.mpr
      (by simpa using Nat.le_of_not_gt hindex)).symm
termination_by factors.size - index
decreasing_by simp_wf; omega

/-- Reading the full source-array interval from zero yields exactly the
array's concrete list, not merely a list with the same product. -/
theorem henselFactorRangeList_full (factors : Array SparsePolyZp) :
    henselFactorRangeList factors factors.size 0 = factors.toList := by
  simpa using henselFactorRangeList_suffix factors 0

theorem henselFactorRangeList_split
    (factors : Array SparsePolyZp) (start mid stop : Nat)
    (hstartMid : start ≤ mid) (hmidStop : mid ≤ stop) :
    henselFactorRangeList factors stop start =
      henselFactorRangeList factors mid start ++
        henselFactorRangeList factors stop mid := by
  rw [henselFactorRangeList]
  by_cases hstartStop : start < stop
  · rw [if_pos hstartStop]
    by_cases hstartMidStrict : start < mid
    · have hmidRange : henselFactorRangeList factors mid start =
          factors[start]! :: henselFactorRangeList factors mid (start + 1) := by
        rw [henselFactorRangeList, if_pos hstartMidStrict]
      rw [hmidRange, List.cons_append]
      congr 1
      exact henselFactorRangeList_split factors (start + 1) mid stop
        (by omega) hmidStop
    · have hstartEqMid : start = mid := by omega
      subst mid
      have hstopRange : henselFactorRangeList factors stop start =
          factors[start]! :: henselFactorRangeList factors stop (start + 1) := by
        rw [henselFactorRangeList, if_pos hstartStop]
      rw [henselFactorRangeList_self, List.nil_append]
      exact hstopRange.symm
  · have hstartEqStop : start = stop := by omega
    subst stop
    have hmidEq : mid = start := by omega
    subst mid
    simp
termination_by stop - start
decreasing_by simp_wf; omega

theorem henselFactorRangeList_singleton
    (factors : Array SparsePolyZp) (start stop : Nat)
    (hlength : stop - start = 1) :
    henselFactorRangeList factors stop start = [factors[start]!] := by
  have hstop : stop = start + 1 := by omega
  subst stop
  rw [henselFactorRangeList, if_pos (by omega)]
  simp

/-- One factor is coprime to a product over a later half-open interval when
it is coprime to every concrete factor in that interval. -/
private theorem henselFactor_coprime_rangeProduct
    (p : Nat) (factors : Array SparsePolyZp)
    (hpairwise : ∀ i j (hi : i < factors.size) (hj : j < factors.size),
      i < j →
      IsCoprime (CLPoly.Math.SparsePolyZp.toPoly p (getElem factors i hi))
        (CLPoly.Math.SparsePolyZp.toPoly p (getElem factors j hj)))
    (source start stop : Nat) (hsource : source < start)
    (hstop : stop ≤ factors.size) :
    IsCoprime (CLPoly.Math.SparsePolyZp.toPoly p factors[source]!)
      (henselFactorRangeProduct p factors stop start) := by
  rw [henselFactorRangeProduct]
  by_cases hmore : start < stop
  · rw [if_pos hmore]
    have hstartSize : start < factors.size := lt_of_lt_of_le hmore hstop
    rw [getElem!_pos factors start hstartSize]
    have hsourceSize : source < factors.size := by omega
    have htail := henselFactor_coprime_rangeProduct p factors hpairwise source
      (start + 1) stop (by omega) hstop
    rw [getElem!_pos factors source hsourceSize] at htail
    rw [getElem!_pos factors source hsourceSize]
    exact (hpairwise source start hsourceSize hstartSize hsource).mul_right htail
  · rw [if_neg hmore]
    exact isCoprime_one_right
termination_by stop - start
decreasing_by simp_wf; omega

/-- Pairwise coprimality of the concrete input array implies coprimality of
the products on every two adjacent half-open intervals. -/
theorem henselFactorRangeProducts_isCoprime
    (p : Nat) (factors : Array SparsePolyZp)
    (hpairwise : ∀ i j (hi : i < factors.size) (hj : j < factors.size),
      i < j →
      IsCoprime (CLPoly.Math.SparsePolyZp.toPoly p (getElem factors i hi))
        (CLPoly.Math.SparsePolyZp.toPoly p (getElem factors j hj)))
    (start mid stop : Nat) (hstartMid : start ≤ mid)
    (hmidStop : mid ≤ stop) (hstop : stop ≤ factors.size) :
    IsCoprime (henselFactorRangeProduct p factors mid start)
      (henselFactorRangeProduct p factors stop mid) := by
  rw [henselFactorRangeProduct]
  by_cases hmore : start < mid
  · rw [if_pos hmore]
    have hstartSize : start < factors.size := by omega
    rw [getElem!_pos factors start hstartSize]
    have hhead := henselFactor_coprime_rangeProduct p factors hpairwise start mid
      stop hmore hstop
    rw [getElem!_pos factors start hstartSize] at hhead
    exact hhead.mul_left
        (henselFactorRangeProducts_isCoprime p factors hpairwise
          (start + 1) mid stop (by omega) hmidStop hstop)
  · rw [if_neg hmore]
    exact isCoprime_one_left
termination_by mid - start
decreasing_by simp_wf; omega

/-- Number of internal nodes created by the source builder on a half-open
factor interval.  Singleton intervals are represented directly by their
parent's `g` or `h`, so only intervals of length at least two allocate nodes. -/
def henselTreeInternalNodeCount : Nat → Nat → Nat
  | start, stop =>
      if hlength : 2 ≤ stop - start then
        let mid := (start + stop) / 2
        1 + henselTreeInternalNodeCount start mid +
          henselTreeInternalNodeCount mid stop
      else 0
termination_by start stop => stop - start
decreasing_by
  · have := Generated.StrictHensel.henselTreeMidpoint_lt_stop start stop hlength
    omega
  · have := Generated.StrictHensel.henselTreeMidpoint_gt_start start stop hlength
    omega

/-- Canonical finite topology induced solely by the source interval split and
preorder append order.  It contains no polynomial values or expected output. -/
def henselTreeBuildTopology : Nat → Nat → Nat →
    Generated.StrictHensel.HenselLiftTree
  | start, stop, root =>
      let mid := (start + stop) / 2
      let leftCount := henselTreeInternalNodeCount start mid
      let left := if 2 ≤ mid - start then
        some (henselTreeBuildTopology start mid (root + 1))
      else none
      let right := if 2 ≤ stop - mid then
        some (henselTreeBuildTopology mid stop (root + 1 + leftCount))
      else none
      .node root left right
termination_by start stop root => stop - start
decreasing_by
  · have hlength : 2 ≤ stop - start := by omega
    have := Generated.StrictHensel.henselTreeMidpoint_lt_stop start stop hlength
    omega
  · have hlength : 2 ≤ stop - start := by omega
    have := Generated.StrictHensel.henselTreeMidpoint_gt_start start stop hlength
    omega

theorem henselTreeInternalNodeCountRawIR_eq
    (start stop : Nat) :
    Generated.StrictHensel.henselTreeInternalNodeCountRawIR start stop =
      henselTreeInternalNodeCount start stop := by
  rw [Generated.StrictHensel.henselTreeInternalNodeCountRawIR,
    henselTreeInternalNodeCount]
  split
  · dsimp only
    rw [henselTreeInternalNodeCountRawIR_eq start ((start + stop) / 2),
      henselTreeInternalNodeCountRawIR_eq ((start + stop) / 2) stop]
  · rfl
termination_by stop - start
decreasing_by all_goals omega

theorem henselTreeBuildTopologyRawIR_eq
    (start stop root : Nat) :
    Generated.StrictHensel.henselTreeBuildTopologyRawIR start stop root =
      henselTreeBuildTopology start stop root := by
  rw [Generated.StrictHensel.henselTreeBuildTopologyRawIR,
    henselTreeBuildTopology]
  rw [henselTreeInternalNodeCountRawIR_eq]
  congr 1
  · split <;> rename_i hbranch
    · congr 1
      exact henselTreeBuildTopologyRawIR_eq start ((start + stop) / 2)
        (root + 1)
    · rfl
  · split <;> rename_i hbranch
    · congr 1
      exact henselTreeBuildTopologyRawIR_eq ((start + stop) / 2) stop
        (root + 1 + henselTreeInternalNodeCount start ((start + stop) / 2))
    · rfl
termination_by stop - start
decreasing_by all_goals omega

@[simp] theorem henselTreeBuildTopology_rootIndex
    (start stop root : Nat) :
    (henselTreeBuildTopology start stop root).rootIndex = root := by
  rw [henselTreeBuildTopology]
  rfl

theorem henselTreeInternalNodeCount_eq
    (start stop : Nat) (hlength : 2 ≤ stop - start) :
    henselTreeInternalNodeCount start stop =
      1 + henselTreeInternalNodeCount start ((start + stop) / 2) +
        henselTreeInternalNodeCount ((start + stop) / 2) stop := by
  rw [henselTreeInternalNodeCount, dif_pos hlength]

theorem henselTreeBuildTopology_nodeCount
    (start stop root : Nat) (hlength : 2 ≤ stop - start) :
    (henselTreeBuildTopology start stop root).nodeCount =
      henselTreeInternalNodeCount start stop := by
  rw [henselTreeInternalNodeCount_eq start stop hlength]
  let mid := (start + stop) / 2
  have hstartMid :=
    Generated.StrictHensel.henselTreeMidpoint_gt_start start stop hlength
  have hmidStop :=
    Generated.StrictHensel.henselTreeMidpoint_lt_stop start stop hlength
  by_cases hleft : 2 ≤ mid - start
  · have hleftCount := henselTreeBuildTopology_nodeCount start mid
      (root + 1) hleft
    by_cases hright : 2 ≤ stop - mid
    · have hrightCount := henselTreeBuildTopology_nodeCount mid stop
        (root + 1 + henselTreeInternalNodeCount start mid) hright
      dsimp [mid] at hleftCount hrightCount
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright ⊢
      unfold Generated.StrictHensel.HenselLiftTree.nodeCount at hleftCount hrightCount ⊢
      simp [hleft, hright, hleftCount, hrightCount]
    · have hrightCount : henselTreeInternalNodeCount mid stop = 0 := by
        rw [henselTreeInternalNodeCount, dif_neg hright]
      dsimp [mid] at hleftCount hrightCount
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright ⊢
      unfold Generated.StrictHensel.HenselLiftTree.nodeCount at hleftCount ⊢
      simp [hleft, hright, hleftCount, hrightCount]
  · have hleftCount : henselTreeInternalNodeCount start mid = 0 := by
      rw [henselTreeInternalNodeCount, dif_neg hleft]
    by_cases hright : 2 ≤ stop - mid
    · have hrightTopologyCount := henselTreeBuildTopology_nodeCount mid stop
        (root + 1 + henselTreeInternalNodeCount start mid) hright
      dsimp [mid] at hleftCount hrightTopologyCount
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright ⊢
      unfold Generated.StrictHensel.HenselLiftTree.nodeCount at hrightTopologyCount ⊢
      rw [hleftCount] at hrightTopologyCount
      simpa [hleft, hright, hleftCount] using hrightTopologyCount
    · have hrightCount : henselTreeInternalNodeCount mid stop = 0 := by
        rw [henselTreeInternalNodeCount, dif_neg hright]
      dsimp [mid] at hleftCount hrightCount
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright ⊢
      unfold Generated.StrictHensel.HenselLiftTree.nodeCount
      simp [hleft, hright, hleftCount, hrightCount]
termination_by stop - start
decreasing_by
  · omega
  · omega
  · omega

/-- The canonical preorder topology uses each mutable-array index exactly
once, and its indices occupy the half-open block beginning at `root`.  This
is the concrete separation fact needed to compose left and right generated
Hensel traversals. -/
theorem henselTreeBuildTopology_indices_nodup_bounded
    (start stop root : Nat) (hlength : 2 ≤ stop - start) :
    let indices := henselLiftTreeIndices
      (henselTreeBuildTopology start stop root)
    indices.Nodup ∧ ∀ index ∈ indices,
      root ≤ index ∧
        index < root + henselTreeInternalNodeCount start stop := by
  let mid := (start + stop) / 2
  have hcount := henselTreeInternalNodeCount_eq start stop hlength
  have hstartMid :=
    Generated.StrictHensel.henselTreeMidpoint_gt_start start stop hlength
  have hmidStop :=
    Generated.StrictHensel.henselTreeMidpoint_lt_stop start stop hlength
  by_cases hleft : 2 ≤ mid - start
  · have hleftIH := henselTreeBuildTopology_indices_nodup_bounded start mid
      (root + 1) hleft
    by_cases hright : 2 ≤ stop - mid
    · have hrightIH := henselTreeBuildTopology_indices_nodup_bounded mid stop
        (root + 1 + henselTreeInternalNodeCount start mid) hright
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright hleftIH hrightIH hcount ⊢
      simp only [hleft, hright, ↓reduceIte, henselLiftTreeIndices]
      constructor
      · apply List.nodup_cons.mpr
        constructor
        · intro hrootMem
          rcases List.mem_append.mp hrootMem with hrootLeft | hrootRight
          · have := hleftIH.2 root hrootLeft
            omega
          · have := hrightIH.2 root hrootRight
            omega
        · exact hleftIH.1.append hrightIH.1 (by
            intro index hindexLeft hindexRight
            have hleftBound := hleftIH.2 index hindexLeft
            have hrightBound := hrightIH.2 index hindexRight
            omega)
      · intro index hindex
        rcases List.mem_cons.mp hindex with rfl | htail
        · omega
        · rcases List.mem_append.mp htail with hindexLeft | hindexRight
          · have := hleftIH.2 index hindexLeft
            omega
          · have := hrightIH.2 index hindexRight
            omega
    · have hrightCount : henselTreeInternalNodeCount mid stop = 0 := by
        rw [henselTreeInternalNodeCount, dif_neg hright]
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright hleftIH hrightCount hcount ⊢
      simp only [hleft, hright, ↓reduceIte, henselLiftTreeIndices]
      constructor
      · rw [List.nodup_cons]
        exact ⟨fun hrootMem => by
          have := hleftIH.2 root hrootMem
          omega, hleftIH.1⟩
      · intro index hindex
        rcases List.mem_cons.mp hindex with rfl | htail
        · omega
        · have := hleftIH.2 index htail
          omega
  · have hleftCount : henselTreeInternalNodeCount start mid = 0 := by
      rw [henselTreeInternalNodeCount, dif_neg hleft]
    by_cases hright : 2 ≤ stop - mid
    · have hrightIH := henselTreeBuildTopology_indices_nodup_bounded mid stop
        (root + 1 + henselTreeInternalNodeCount start mid) hright
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright hrightIH hleftCount hcount ⊢
      simp only [hleft, hright, ↓reduceIte, henselLiftTreeIndices]
      constructor
      · rw [List.nodup_cons]
        exact ⟨fun hrootMem => by
          have := hrightIH.2 root hrootMem
          omega, hrightIH.1⟩
      · intro index hindex
        rcases List.mem_cons.mp hindex with rfl | htail
        · omega
        · have := hrightIH.2 index htail
          omega
    · have hrightCount : henselTreeInternalNodeCount mid stop = 0 := by
        rw [henselTreeInternalNodeCount, dif_neg hright]
      rw [henselTreeBuildTopology]
      dsimp [mid] at hleft hright hleftCount hrightCount hcount ⊢
      simp only [hleft, hright, ↓reduceIte, henselLiftTreeIndices]
      constructor
      · simp
      · intro index hindex
        have hindexEq : index = root := by simpa using hindex
        rw [hindexEq]
        omega
termination_by stop - start
decreasing_by all_goals omega

theorem nat_toUInt32_toInt32_toInt_eq (n : Nat) (hbound : n < 2 ^ 31) :
    n.toUInt32.toInt32.toInt = (n : Int) := by
  change (Int32.ofNat n).toInt = (n : Int)
  exact Int32.toInt_ofNat_of_lt hbound

/-- Exact raw pointer values installed at one node by the source preorder
builder, before interpreting their signed integer representation. -/
def HenselTreeNodeRawTopologyMatches (start stop parent : Nat)
    (node : HenselNode) : Prop :=
  let mid := (start + stop) / 2
  node.left =
      (if 2 ≤ mid - start then (parent + 1).toUInt32.toInt32 else -1) ∧
  node.right =
      (if 2 ≤ stop - mid then
        (parent + 1 + henselTreeInternalNodeCount start mid).toUInt32.toInt32
      else -1)

theorem HenselTreeNodeRawTopologyMatches.liftChildMatches
    {start stop parent : Nat} {node : HenselNode}
    (hlength : 2 ≤ stop - start)
    (hraw : HenselTreeNodeRawTopologyMatches start stop parent node)
    (hbound : parent + henselTreeInternalNodeCount start stop < 2 ^ 31) :
    let tree := henselTreeBuildTopology start stop parent
    liftChildMatches node.left (match tree with | .node _ left _ => left) ∧
    liftChildMatches node.right (match tree with | .node _ _ right => right) := by
  let mid := (start + stop) / 2
  have hcount := henselTreeInternalNodeCount_eq start stop hlength
  simp only [HenselTreeNodeRawTopologyMatches] at hraw
  rw [henselTreeBuildTopology]
  dsimp only
  by_cases hleft : 2 ≤ mid - start
  · have hleftPositive : 0 < henselTreeInternalNodeCount start mid := by
      rw [henselTreeInternalNodeCount_eq start mid hleft]
      omega
    by_cases hright : 2 ≤ stop - mid
    · have hrightPositive : 0 < henselTreeInternalNodeCount mid stop := by
        rw [henselTreeInternalNodeCount_eq mid stop hright]
        omega
      dsimp [mid] at hleft hright hraw ⊢
      rw [if_pos hleft, if_pos hright] at hraw ⊢
      constructor
      · rw [hraw.1]
        change (parent + 1).toUInt32.toInt32.toInt = _
        rw [henselTreeBuildTopology_rootIndex,
          nat_toUInt32_toInt32_toInt_eq]
        omega
      · rw [hraw.2]
        change (parent + 1 + henselTreeInternalNodeCount start
          ((start + stop) / 2)).toUInt32.toInt32.toInt = _
        rw [henselTreeBuildTopology_rootIndex,
          nat_toUInt32_toInt32_toInt_eq]
        dsimp [mid] at hcount hbound ⊢
        omega
    · dsimp [mid] at hleft hright hraw ⊢
      rw [if_pos hleft, if_neg hright] at hraw ⊢
      constructor
      · rw [hraw.1]
        change (parent + 1).toUInt32.toInt32.toInt = _
        rw [henselTreeBuildTopology_rootIndex,
          nat_toUInt32_toInt32_toInt_eq]
        dsimp [mid] at hcount hbound ⊢
        omega
      · exact hraw.2
  · by_cases hright : 2 ≤ stop - mid
    · have hrightPositive : 0 < henselTreeInternalNodeCount mid stop := by
        rw [henselTreeInternalNodeCount_eq mid stop hright]
        omega
      dsimp [mid] at hleft hright hraw ⊢
      rw [if_neg hleft, if_pos hright] at hraw ⊢
      constructor
      · exact hraw.1
      · rw [hraw.2]
        change (parent + 1 + henselTreeInternalNodeCount start
          ((start + stop) / 2)).toUInt32.toInt32.toInt = _
        rw [henselTreeBuildTopology_rootIndex,
          nat_toUInt32_toInt32_toInt_eq]
        dsimp [mid] at hcount hbound ⊢
        omega
    · dsimp [mid] at hleft hright hraw ⊢
      rw [if_neg hleft, if_neg hright] at hraw ⊢
      exact hraw

/-- Mathematical content stored at one freshly constructed tree node before
pairwise coprimality specializes its gcd to one. -/
noncomputable def HenselTreeNodeGCDInvariant (p : Nat) [Fact (Nat.Prime p)]
    (factors : Array SparsePolyZp) (start stop : Nat)
    (node : HenselNode) : Prop :=
  let mid := (start + stop) / 2
  let left := henselFactorRangeProduct p factors mid start
  let right := henselFactorRangeProduct p factors stop mid
  toPolyMod p node.g = left ∧
  toPolyMod p node.h = right ∧
  toPolyMod p node.s * left + toPolyMod p node.t * right =
      GCDMonoid.gcd left right ∧
  (GCDMonoid.gcd left right).Monic

/-- Initial Hensel invariant at a tree node.  Unlike the gcd form above this
is the exact certificate consumed by the later lifting pass. -/
noncomputable def HenselTreeNodeInitialInvariant (p : Nat)
    (factors : Array SparsePolyZp) (start stop : Nat)
    (node : HenselNode) : Prop :=
  let mid := (start + stop) / 2
  let left := henselFactorRangeProduct p factors mid start
  let right := henselFactorRangeProduct p factors stop mid
  toPolyMod p node.g = left ∧
  toPolyMod p node.h = right ∧
  toPolyMod p node.s * left + toPolyMod p node.t * right = 1

theorem sparsePolyZpToPoly_ne_zero_of_canonical_nonempty
    (p : Nat) [Fact (Nat.Prime p)] (f : SparsePolyZp)
    (hcanonical : CLPoly.Math.SparsePolyZp.Canonical p f)
    (hnonempty : 0 < f.size) :
    CLPoly.Math.SparsePolyZp.toPoly p f ≠ 0 := by
  intro hzero
  have hdegree :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      p f hcanonical hnonempty
  rw [hzero, Polynomial.degree_zero] at hdegree
  simp at hdegree

private theorem strictHenselTreeProductLoopRawIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (stop : Nat) (hstop : stop ≤ factors.size) :
    ∀ index product,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat product →
      0 < product.size →
      ∃ output,
        Generated.StrictHensel.henselTreeProductLoopRawIR
            (strictHenselTreeBuildRawOps this mulProvider) factors stop index
            product = .ok output ∧
        CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output ∧
        0 < output.size ∧
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output =
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat product *
            henselFactorRangeProduct this._p.toNat factors stop index := by
  intro index product hproduct hproductNonempty
  rw [Generated.StrictHensel.henselTreeProductLoopRawIR]
  rw [henselFactorRangeProduct]
  by_cases hmore : index < stop
  · rw [dif_pos hmore]
    have hindex : index < factors.size := lt_of_lt_of_le hmore hstop
    rw [getElem!_pos factors index hindex]
    have hget : factors[index]? = some factors[index] :=
      Array.getElem?_eq_getElem hindex
    rw [hget]
    have hfactorMem : factors[index] ∈ factors.toList :=
      Array.getElem_mem_toList hindex
    rcases StrictDDF.strictMulIR_refines_mul this hcfg product factors[index]
        mulProvider hproduct (hfactors factors[index] hfactorMem) with
      ⟨next, hnextRun, hnextCanonical, hnextPoly⟩
    have hproductPolyNonzero :=
      sparsePolyZpToPoly_ne_zero_of_canonical_nonempty this._p.toNat product
        hproduct hproductNonempty
    have hfactorPolyNonzero :=
      sparsePolyZpToPoly_ne_zero_of_canonical_nonempty this._p.toNat
        factors[index] (hfactors factors[index] hfactorMem)
        (hfactorsNonempty factors[index] hfactorMem)
    have hnextPolyNonzero :
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat next ≠ 0 := by
      rw [hnextPoly]
      exact mul_ne_zero hproductPolyNonzero hfactorPolyNonzero
    have hnextNonempty : 0 < next.size := by
      by_contra hnot
      have hsize : next.size = 0 := Nat.eq_zero_of_not_pos hnot
      have hempty : next = #[] := Array.size_eq_zero_iff.mp hsize
      subst next
      exact hnextPolyNonzero
        (CLPoly.Math.SparsePolyZp.toPoly_empty this._p.toNat)
    simp only [strictHenselTreeBuildRawOps, hnextRun, bind, Except.bind]
    rcases strictHenselTreeProductLoopRawIR_refines this hcfg mulProvider
        factors hfactors hfactorsNonempty stop hstop (index + 1) next
        hnextCanonical hnextNonempty with
      ⟨output, houtputRun, houtputCanonical, houtputNonempty, houtputPoly⟩
    refine ⟨output, houtputRun, houtputCanonical, houtputNonempty, ?_⟩
    rw [houtputPoly, hnextPoly]
    rw [if_pos hmore]
    ring
  · rw [dif_neg hmore]
    exact ⟨product, rfl, hproduct, hproductNonempty, by
      rw [if_neg hmore]
      ring⟩
termination_by index product _ _ => stop - index
decreasing_by simp_wf; omega

/-- Total semantic contract for each of the two actual factor-product loops
used by the Hensel tree builder. -/
theorem strictHenselTreeProductRangeRawIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (start stop : Nat) (hnonempty : start < stop)
    (hstop : stop ≤ factors.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeProductRangeRawIR
          (strictHenselTreeBuildRawOps this mulProvider) factors start stop =
        .ok output ∧
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output ∧
      0 < output.size ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output =
        henselFactorRangeProduct this._p.toNat factors stop start := by
  have hstart : start < factors.size := lt_of_lt_of_le hnonempty hstop
  have hget : factors[start]? = some factors[start] :=
    Array.getElem?_eq_getElem hstart
  have hfirstCanonical := hfactors factors[start]
    (Array.getElem_mem_toList hstart)
  rcases strictHenselTreeProductLoopRawIR_refines this hcfg mulProvider
      factors hfactors hfactorsNonempty stop hstop (start + 1)
      factors[start] hfirstCanonical
      (hfactorsNonempty factors[start] (Array.getElem_mem_toList hstart)) with
    ⟨output, hrun, hcanonical, houtputNonempty, hpoly⟩
  refine ⟨output, ?_, hcanonical, houtputNonempty, ?_⟩
  · simp [Generated.StrictHensel.henselTreeProductRangeRawIR, hnonempty,
      hget, hrun]
  · rw [henselFactorRangeProduct, if_pos hnonempty, hpoly]
    rw [getElem!_pos factors start hstart]

theorem henselTreeSetNodeRawIR_succeeds
    (nodes : Array HenselNode) (index : Nat)
    (update : HenselNode → HenselNode) (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeSetNodeRawIR nodes index update =
        .ok output ∧
      output.size = nodes.size := by
  refine ⟨nodes.set! index (update nodes[index]), ?_, by simp⟩
  simp [Generated.StrictHensel.henselTreeSetNodeRawIR,
    Array.getElem?_eq_getElem hindex]

/-- A checked node update has the ordinary array-frame behavior: it changes
the selected element exactly as requested and leaves every other element
unchanged. -/
theorem henselTreeSetNodeRawIR_frame
    (nodes : Array HenselNode) (index : Nat)
    (update : HenselNode → HenselNode) (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeSetNodeRawIR nodes index update =
        .ok output ∧
      output.size = nodes.size ∧
      ∃ houtput : index < output.size,
      getElem output index houtput = update nodes[index] ∧
      ∀ other, ∀ hother : other < nodes.size, other ≠ index →
        ∀ hotherOutput : other < output.size,
        getElem output other hotherOutput = getElem nodes other hother := by
  let output := nodes.set! index (update nodes[index])
  have hsize : output.size = nodes.size := by simp [output]
  have houtput : index < output.size := by omega
  refine ⟨output, ?_, hsize, houtput, ?_, ?_⟩
  · simp [Generated.StrictHensel.henselTreeSetNodeRawIR,
      Array.getElem?_eq_getElem hindex, output]
  · simp [output]
  · intro other hother hne hotherOutput
    have hset := Array.getElem_setIfInBounds_ne
      (xs := nodes) (i := index) (j := other)
      (a := update nodes[index]) hother hne.symm
    simpa [output] using hset

theorem henselTreeStoreNodeRawIR_succeeds
    (nodes : Array HenselNode) (index : Nat)
    (g h s t : SparsePolyZp) (start stop : Nat)
    (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      output.size = nodes.size := by
  rcases henselTreeSetNodeRawIR_succeeds nodes index
      (fun node => { node with g :=
        Generated.StrictHensel.henselTreeZpToZZIR g }) hindex with
    ⟨nodes1, hnodes1, hsize1⟩
  have hindex1 : index < nodes1.size := by omega
  rcases henselTreeSetNodeRawIR_succeeds nodes1 index
      (fun node => { node with h :=
        Generated.StrictHensel.henselTreeZpToZZIR h }) hindex1 with
    ⟨nodes2, hnodes2, hsize2⟩
  have hindex2 : index < nodes2.size := by omega
  rcases henselTreeSetNodeRawIR_succeeds nodes2 index
      (fun node => { node with s :=
        Generated.StrictHensel.henselTreeZpToZZIR s }) hindex2 with
    ⟨nodes3, hnodes3, hsize3⟩
  have hindex3 : index < nodes3.size := by omega
  rcases henselTreeSetNodeRawIR_succeeds nodes3 index
      (fun node => { node with t :=
        Generated.StrictHensel.henselTreeZpToZZIR t }) hindex3 with
    ⟨nodes4, hnodes4, hsize4⟩
  have hindex4 : index < nodes4.size := by omega
  rcases henselTreeSetNodeRawIR_succeeds nodes4 index
      (fun node => { node with leaf_start := start.toUInt32.toInt32 })
      hindex4 with ⟨nodes5, hnodes5, hsize5⟩
  have hindex5 : index < nodes5.size := by omega
  rcases henselTreeSetNodeRawIR_succeeds nodes5 index
      (fun node => { node with leaf_end := stop.toUInt32.toInt32 })
      hindex5 with ⟨output, houtput, hsize6⟩
  refine ⟨output, ?_, by omega⟩
  simp only [Generated.StrictHensel.henselTreeStoreNodeRawIR, hnodes1,
    hnodes2, hnodes3, hnodes4, hnodes5, houtput, bind, Except.bind]

/-- The six stores affect no array entry other than the selected node. -/
theorem henselTreeStoreNodeRawIR_frame
    (nodes : Array HenselNode) (index : Nat)
    (g h s t : SparsePolyZp) (start stop : Nat)
    (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      output.size = nodes.size ∧
      ∀ other, ∀ hother : other < nodes.size, other ≠ index →
        ∀ hotherOutput : other < output.size,
        getElem output other hotherOutput = getElem nodes other hother := by
  rcases henselTreeSetNodeRawIR_frame nodes index
      (fun node => { node with g :=
        Generated.StrictHensel.henselTreeZpToZZIR g }) hindex with
    ⟨nodes1, hrun1, hsize1, hindex1, _, hframe1⟩
  rcases henselTreeSetNodeRawIR_frame nodes1 index
      (fun node => { node with h :=
        Generated.StrictHensel.henselTreeZpToZZIR h }) hindex1 with
    ⟨nodes2, hrun2, hsize2, hindex2, _, hframe2⟩
  rcases henselTreeSetNodeRawIR_frame nodes2 index
      (fun node => { node with s :=
        Generated.StrictHensel.henselTreeZpToZZIR s }) hindex2 with
    ⟨nodes3, hrun3, hsize3, hindex3, _, hframe3⟩
  rcases henselTreeSetNodeRawIR_frame nodes3 index
      (fun node => { node with t :=
        Generated.StrictHensel.henselTreeZpToZZIR t }) hindex3 with
    ⟨nodes4, hrun4, hsize4, hindex4, _, hframe4⟩
  rcases henselTreeSetNodeRawIR_frame nodes4 index
      (fun node => { node with leaf_start := start.toUInt32.toInt32 })
      hindex4 with
    ⟨nodes5, hrun5, hsize5, hindex5, _, hframe5⟩
  rcases henselTreeSetNodeRawIR_frame nodes5 index
      (fun node => { node with leaf_end := stop.toUInt32.toInt32 })
      hindex5 with
    ⟨output, hrun6, hsize6, _, _, hframe6⟩
  refine ⟨output, ?_, by omega, ?_⟩
  · simp only [Generated.StrictHensel.henselTreeStoreNodeRawIR, hrun1,
      hrun2, hrun3, hrun4, hrun5, hrun6, bind, Except.bind]
  intro other hother hne hotherOutput
  have hother1 : other < nodes1.size := by omega
  have hother2 : other < nodes2.size := by omega
  have hother3 : other < nodes3.size := by omega
  have hother4 : other < nodes4.size := by omega
  have hother5 : other < nodes5.size := by omega
  calc
    getElem output other hotherOutput = getElem nodes5 other hother5 :=
      hframe6 other hother5 hne hotherOutput
    _ = getElem nodes4 other hother4 := hframe5 other hother4 hne hother5
    _ = getElem nodes3 other hother3 := hframe4 other hother3 hne hother4
    _ = getElem nodes2 other hother2 := hframe3 other hother2 hne hother3
    _ = getElem nodes1 other hother1 := hframe2 other hother1 hne hother2
    _ = getElem nodes other hother := hframe1 other hother hne hother1

/-- Exact semantic effect of the six checked writes at a tree node. -/
theorem henselTreeStoreNodeRawIR_refines
    (nodes : Array HenselNode) (index : Nat)
    (g h s t : SparsePolyZp) (start stop : Nat)
    (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      output.size = nodes.size ∧
      ∃ houtput : index < output.size,
      (getElem output index houtput).g =
          Generated.StrictHensel.henselTreeZpToZZIR g ∧
      (getElem output index houtput).h =
          Generated.StrictHensel.henselTreeZpToZZIR h ∧
      (getElem output index houtput).s =
          Generated.StrictHensel.henselTreeZpToZZIR s ∧
      (getElem output index houtput).t =
          Generated.StrictHensel.henselTreeZpToZZIR t ∧
      (getElem output index houtput).leaf_start = start.toUInt32.toInt32 ∧
      (getElem output index houtput).leaf_end = stop.toUInt32.toInt32 := by
  rcases henselTreeSetNodeRawIR_frame nodes index
      (fun node => { node with g :=
        Generated.StrictHensel.henselTreeZpToZZIR g }) hindex with
    ⟨nodes1, hrun1, hsize1, hindex1, hselected1, _⟩
  rcases henselTreeSetNodeRawIR_frame nodes1 index
      (fun node => { node with h :=
        Generated.StrictHensel.henselTreeZpToZZIR h }) hindex1 with
    ⟨nodes2, hrun2, hsize2, hindex2, hselected2, _⟩
  rcases henselTreeSetNodeRawIR_frame nodes2 index
      (fun node => { node with s :=
        Generated.StrictHensel.henselTreeZpToZZIR s }) hindex2 with
    ⟨nodes3, hrun3, hsize3, hindex3, hselected3, _⟩
  rcases henselTreeSetNodeRawIR_frame nodes3 index
      (fun node => { node with t :=
        Generated.StrictHensel.henselTreeZpToZZIR t }) hindex3 with
    ⟨nodes4, hrun4, hsize4, hindex4, hselected4, _⟩
  rcases henselTreeSetNodeRawIR_frame nodes4 index
      (fun node => { node with leaf_start := start.toUInt32.toInt32 })
      hindex4 with
    ⟨nodes5, hrun5, hsize5, hindex5, hselected5, _⟩
  rcases henselTreeSetNodeRawIR_frame nodes5 index
      (fun node => { node with leaf_end := stop.toUInt32.toInt32 })
      hindex5 with
    ⟨output, hrun6, hsize6, houtput, hselected6, _⟩
  refine ⟨output, ?_, by omega, houtput, ?_⟩
  · simp only [Generated.StrictHensel.henselTreeStoreNodeRawIR, hrun1,
      hrun2, hrun3, hrun4, hrun5, hrun6, bind, Except.bind]
  simp [hselected6, hselected5, hselected4, hselected3, hselected2,
    hselected1]

/-- The six checked writes of the concrete tree builder store canonical
integer `g`/`h` fields whenever their generated finite-field sources are
canonical. -/
theorem henselTreeStoreNodeRawIR_canonical
    (p : Nat) (nodes : Array HenselNode) (index : Nat)
    (g h s t : SparsePolyZp) (start stop : Nat)
    (hindex : index < nodes.size)
    (hgCanonical : CLPoly.Math.SparsePolyZp.Canonical p g)
    (hhCanonical : CLPoly.Math.SparsePolyZp.Canonical p h) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      output.size = nodes.size ∧
      ∃ houtput : index < output.size,
        HenselNodeCanonical (getElem output index houtput) := by
  rcases henselTreeStoreNodeRawIR_refines nodes index g h s t start stop
      hindex with
    ⟨output, hrun, hsize, houtput, houtputG, houtputH, _⟩
  refine ⟨output, hrun, hsize, houtput, ?_⟩
  unfold HenselNodeCanonical
  rw [houtputG, houtputH]
  exact ⟨henselTreeZpToZZIR_canonical p g hgCanonical,
    henselTreeZpToZZIR_canonical p h hhCanonical⟩

/-- The actual six-write builder store exposes both canonicality and the
literal coefficient-one head of its converted right interval product. -/
theorem henselTreeStoreNodeRawIR_canonical_oneHead
    (p : Nat) [Fact (Nat.Prime p)] (nodes : Array HenselNode) (index : Nat)
    (g h s t : SparsePolyZp) (start stop : Nat)
    (hindex : index < nodes.size)
    (hgCanonical : CLPoly.Math.SparsePolyZp.Canonical p g)
    (hhCanonical : CLPoly.Math.SparsePolyZp.Canonical p h)
    (hhNonempty : 0 < h.size)
    (hhMonic : (CLPoly.Math.SparsePolyZp.toPoly p h).Monic) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      output.size = nodes.size ∧
      ∃ houtput : index < output.size,
        HenselNodeCanonical (getElem output index houtput) ∧
        HasPhysicalOneHead (getElem output index houtput).h := by
  rcases henselTreeStoreNodeRawIR_refines nodes index g h s t start stop
      hindex with
    ⟨output, hrun, hsize, houtput, houtputG, houtputH, _⟩
  refine ⟨output, hrun, hsize, houtput, ?_, ?_⟩
  · unfold HenselNodeCanonical
    rw [houtputG, houtputH]
    exact ⟨henselTreeZpToZZIR_canonical p g hgCanonical,
      henselTreeZpToZZIR_canonical p h hhCanonical⟩
  · rw [houtputH]
    exact henselTreeZpToZZIR_hasPhysicalOneHead_of_monic p h hhCanonical
      hhNonempty hhMonic

/-- Storing the two concrete range products and the concrete EEA result
establishes the node's exact gcd semantics. -/
theorem henselTreeStoreNodeRawIR_refines_gcd
    (p : Nat) [Fact (Nat.Prime p)] (factors : Array SparsePolyZp)
    (nodes : Array HenselNode)
    (index start stop : Nat) (g h gcd s t : SparsePolyZp)
    (hindex : index < nodes.size)
    (hg : CLPoly.Math.SparsePolyZp.toPoly p g =
      henselFactorRangeProduct p factors ((start + stop) / 2) start)
    (hh : CLPoly.Math.SparsePolyZp.toPoly p h =
      henselFactorRangeProduct p factors stop ((start + stop) / 2))
    (hbezout : CLPoly.Math.SparsePolyZp.toPoly p gcd =
      CLPoly.Math.SparsePolyZp.toPoly p s *
          CLPoly.Math.SparsePolyZp.toPoly p g +
        CLPoly.Math.SparsePolyZp.toPoly p t *
          CLPoly.Math.SparsePolyZp.toPoly p h)
    (hmonic : (CLPoly.Math.SparsePolyZp.toPoly p gcd).Monic)
    (hgcd : CLPoly.Math.SparsePolyZp.toPoly p gcd =
      GCDMonoid.gcd (CLPoly.Math.SparsePolyZp.toPoly p g)
        (CLPoly.Math.SparsePolyZp.toPoly p h)) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      ∃ houtput : index < output.size,
      output.size = nodes.size ∧
      HenselTreeNodeGCDInvariant p factors start stop
        (getElem output index houtput) := by
  rcases henselTreeStoreNodeRawIR_refines nodes index g h s t start stop
      hindex with
    ⟨output, hrun, hsize, houtput, houtputG, houtputH, houtputS, houtputT,
      houtputStart, houtputStop⟩
  refine ⟨output, hrun, houtput, hsize, ?_⟩
  simp only [HenselTreeNodeGCDInvariant]
  rw [houtputG, houtputH, houtputS, houtputT,
    henselTreeZpToZZIR_toPolyMod, henselTreeZpToZZIR_toPolyMod,
    henselTreeZpToZZIR_toPolyMod, henselTreeZpToZZIR_toPolyMod]
  refine ⟨hg, hh, ?_, ?_⟩
  · rw [← hg, ← hh, ← hgcd, hbezout]
  · rw [← hg, ← hh, ← hgcd]
    exact hmonic

/-- When the two interval products are coprime, the same concrete EEA write
establishes the exact unit Bézout certificate required by Hensel lifting. -/
theorem henselTreeStoreNodeRawIR_refines_initial
    (p : Nat) [Fact (Nat.Prime p)] (factors : Array SparsePolyZp)
    (nodes : Array HenselNode)
    (index start stop : Nat) (g h gcd s t : SparsePolyZp)
    (hindex : index < nodes.size)
    (hg : CLPoly.Math.SparsePolyZp.toPoly p g =
      henselFactorRangeProduct p factors ((start + stop) / 2) start)
    (hh : CLPoly.Math.SparsePolyZp.toPoly p h =
      henselFactorRangeProduct p factors stop ((start + stop) / 2))
    (hbezout : CLPoly.Math.SparsePolyZp.toPoly p gcd =
      CLPoly.Math.SparsePolyZp.toPoly p s *
          CLPoly.Math.SparsePolyZp.toPoly p g +
        CLPoly.Math.SparsePolyZp.toPoly p t *
          CLPoly.Math.SparsePolyZp.toPoly p h)
    (hmonic : (CLPoly.Math.SparsePolyZp.toPoly p gcd).Monic)
    (hgcd : CLPoly.Math.SparsePolyZp.toPoly p gcd =
      GCDMonoid.gcd (CLPoly.Math.SparsePolyZp.toPoly p g)
        (CLPoly.Math.SparsePolyZp.toPoly p h))
    (hcoprime : IsCoprime
      (henselFactorRangeProduct p factors ((start + stop) / 2) start)
      (henselFactorRangeProduct p factors stop ((start + stop) / 2))) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      ∃ houtput : index < output.size,
      output.size = nodes.size ∧
      HenselTreeNodeInitialInvariant p factors start stop
        (getElem output index houtput) := by
  rcases henselTreeStoreNodeRawIR_refines_gcd p factors nodes index start stop
      g h gcd s t hindex hg hh hbezout hmonic hgcd with
    ⟨output, hrun, houtput, hsize, hinvariant⟩
  have hcoprime' : IsCoprime
      (CLPoly.Math.SparsePolyZp.toPoly p g)
      (CLPoly.Math.SparsePolyZp.toPoly p h) := by
    rw [hg, hh]
    exact hcoprime
  have hunitGCD : IsUnit
      (GCDMonoid.gcd (CLPoly.Math.SparsePolyZp.toPoly p g)
        (CLPoly.Math.SparsePolyZp.toPoly p h)) :=
    (gcd_isUnit_iff _ _).2 hcoprime'
  have hunit : IsUnit (CLPoly.Math.SparsePolyZp.toPoly p gcd) := by
    rw [hgcd]
    exact hunitGCD
  have honeRaw : CLPoly.Math.SparsePolyZp.toPoly p gcd = 1 :=
    hmonic.eq_one_of_isUnit hunit
  have hone : GCDMonoid.gcd (CLPoly.Math.SparsePolyZp.toPoly p g)
      (CLPoly.Math.SparsePolyZp.toPoly p h) = 1 := by
    rw [← hgcd, honeRaw]
  refine ⟨output, hrun, houtput, hsize, ?_⟩
  refine ⟨hinvariant.1, hinvariant.2.1, ?_⟩
  rw [← hone, hg, hh]
  exact hinvariant.2.2.1

/-- Array frame used by recursive construction: the array may grow and the
current node may change, but every other pre-existing node is identical. -/
def HenselTreeArrayFrame (before after : Array HenselNode)
    (current : Nat) : Prop :=
  before.size ≤ after.size ∧
  ∀ index, index ≠ current →
    ∀ (hbefore : index < before.size) (hafter : index < after.size),
      getElem after index hafter = getElem before index hbefore

/-- Equality of the algebraic payload of a node; child pointers are excluded
because installing them must not invalidate polynomial semantics. -/
def HenselTreeNodeAlgebraEq (left right : HenselNode) : Prop :=
  left.g = right.g ∧ left.h = right.h ∧ left.s = right.s ∧ left.t = right.t

theorem HenselTreeNodeAlgebraEq.refl (node : HenselNode) :
    HenselTreeNodeAlgebraEq node node := by
  exact ⟨rfl, rfl, rfl, rfl⟩

theorem HenselTreeNodeAlgebraEq.of_eq {left right : HenselNode}
    (heq : left = right) : HenselTreeNodeAlgebraEq left right := by
  subst right
  exact HenselTreeNodeAlgebraEq.refl left

theorem HenselTreeNodeAlgebraEq.trans {first second third : HenselNode}
    (hfirst : HenselTreeNodeAlgebraEq first second)
    (hsecond : HenselTreeNodeAlgebraEq second third) :
    HenselTreeNodeAlgebraEq first third := by
  exact ⟨hfirst.1.trans hsecond.1,
    hfirst.2.1.trans hsecond.2.1,
    hfirst.2.2.1.trans hsecond.2.2.1,
    hfirst.2.2.2.trans hsecond.2.2.2⟩

theorem HenselNodeCanonical.of_algebraEq
    {source target : HenselNode} (hsource : HenselNodeCanonical source)
    (heq : HenselTreeNodeAlgebraEq target source) :
    HenselNodeCanonical target := by
  rcases heq with ⟨hg, hh, _, _⟩
  simpa [HenselNodeCanonical, hg, hh] using hsource

theorem HasPhysicalOneHead.of_algebraEq
    {source target : HenselNode} (hsource : HasPhysicalOneHead source.h)
    (heq : HenselTreeNodeAlgebraEq target source) :
    HasPhysicalOneHead target.h := by
  simpa [heq.2.1] using hsource

theorem HenselTreeArrayFrame.canonical_of_size_eq
    {before after : Array HenselNode} {current : Nat}
    (hframe : HenselTreeArrayFrame before after current)
    (hsize : after.size = before.size)
    (hbefore : HenselArrayCanonical before)
    (hcurrent : current < after.size)
    (hvalue : HenselNodeCanonical after[current]) :
    HenselArrayCanonical after := by
  intro index hindex
  by_cases heq : index = current
  · subst index
    exact hvalue
  · have hindexBefore : index < before.size := by omega
    have hsame := hframe.2 index heq hindexBefore hindex
    rw [hsame]
    exact hbefore index hindexBefore

theorem HenselArrayCanonical.push_default
    {nodes : Array HenselNode} (hnodes : HenselArrayCanonical nodes) :
    HenselArrayCanonical (nodes.push default) := by
  intro index hindex
  by_cases hbefore : index < nodes.size
  · simpa only [Array.getElem_push, hbefore, ↓reduceDIte] using
      hnodes index hbefore
  · have hindex' : index < nodes.size + 1 := by simpa using hindex
    have heq : index = nodes.size := by omega
    subst index
    simp only [Array.getElem_push, lt_self_iff_false, ↓reduceDIte]
    change StrictPolynomialMod.SparsePolyZZCanonical (#[] : SparsePolyZZ) ∧
      StrictPolynomialMod.SparsePolyZZCanonical (#[] : SparsePolyZZ)
    simp [StrictPolynomialMod.SparsePolyZZCanonical]

theorem HenselArrayCanonical.singleton_default :
    HenselArrayCanonical #[default] := by
  intro index hindex
  have heq : index = 0 := by simpa using hindex
  subst index
  change StrictPolynomialMod.SparsePolyZZCanonical (#[] : SparsePolyZZ) ∧
    StrictPolynomialMod.SparsePolyZZCanonical (#[] : SparsePolyZZ)
  simp [StrictPolynomialMod.SparsePolyZZCanonical]

theorem HenselTreeNodeGCDInvariant.of_algebraEq
    {p : Nat} [Fact (Nat.Prime p)] {factors : Array SparsePolyZp}
    {start stop : Nat}
    {source target : HenselNode}
    (hinvariant : HenselTreeNodeGCDInvariant p factors start stop source)
    (heq : HenselTreeNodeAlgebraEq target source) :
    HenselTreeNodeGCDInvariant p factors start stop target := by
  rcases heq with ⟨hg, hh, hs, ht⟩
  simpa [HenselTreeNodeGCDInvariant, hg, hh, hs, ht] using hinvariant

theorem HenselTreeNodeGCDInvariant.toInitial
    {p : Nat} [Fact (Nat.Prime p)] {factors : Array SparsePolyZp}
    {start stop : Nat}
    {node : HenselNode}
    (hinvariant : HenselTreeNodeGCDInvariant p factors start stop node)
    (hcoprime : IsCoprime
      (henselFactorRangeProduct p factors ((start + stop) / 2) start)
      (henselFactorRangeProduct p factors stop ((start + stop) / 2))) :
    HenselTreeNodeInitialInvariant p factors start stop node := by
  simp only [HenselTreeNodeGCDInvariant] at hinvariant
  simp only [HenselTreeNodeInitialInvariant]
  have hunit : IsUnit (GCDMonoid.gcd
      (henselFactorRangeProduct p factors ((start + stop) / 2) start)
      (henselFactorRangeProduct p factors stop ((start + stop) / 2))) :=
    (gcd_isUnit_iff _ _).2 hcoprime
  have hone := hinvariant.2.2.2.eq_one_of_isUnit hunit
  exact ⟨hinvariant.1, hinvariant.2.1, by
    simpa [hone] using hinvariant.2.2.1⟩

theorem HenselTreeArrayFrame.refl (nodes : Array HenselNode)
    (current : Nat) : HenselTreeArrayFrame nodes nodes current := by
  exact ⟨by simp, fun _ _ _ _ => rfl⟩

theorem HenselTreeArrayFrame.trans {first second third : Array HenselNode}
    {current : Nat} (hfirst : HenselTreeArrayFrame first second current)
    (hsecond : HenselTreeArrayFrame second third current) :
    HenselTreeArrayFrame first third current := by
  refine ⟨le_trans hfirst.1 hsecond.1, ?_⟩
  intro index hne hfirstIndex hthirdIndex
  have hsecondIndex : index < second.size :=
    lt_of_lt_of_le hfirstIndex hfirst.1
  exact (hsecond.2 index hne hsecondIndex hthirdIndex).trans
    (hfirst.2 index hne hfirstIndex hsecondIndex)

/-- A recursive call at a freshly appended child composes with the parent's
frame because no old index can equal that child index. -/
theorem HenselTreeArrayFrame.trans_fresh
    {first second third : Array HenselNode} {current child : Nat}
    (hfirst : HenselTreeArrayFrame first second current)
    (hsecond : HenselTreeArrayFrame second third child)
    (hfresh : first.size ≤ child) :
    HenselTreeArrayFrame first third current := by
  refine ⟨le_trans hfirst.1 hsecond.1, ?_⟩
  intro index hne hfirstIndex hthirdIndex
  have hchild : index ≠ child := by omega
  have hsecondIndex : index < second.size :=
    lt_of_lt_of_le hfirstIndex hfirst.1
  exact (hsecond.2 index hchild hsecondIndex hthirdIndex).trans
    (hfirst.2 index hne hfirstIndex hsecondIndex)

theorem HenselTreeArrayFrame.preservesFrom
    {before after : Array HenselNode} {current lower : Nat}
    (hframe : HenselTreeArrayFrame before after current)
    (hcurrent : current < lower) :
    HenselTreePreservesFrom lower before after := by
  refine ⟨hframe.1, ?_⟩
  intro index hlower hbefore
  have hafter : index < after.size := lt_of_lt_of_le hbefore hframe.1
  rw [Array.getElem?_eq_getElem hafter,
    Array.getElem?_eq_getElem hbefore]
  exact congrArg some (hframe.2 index (by omega) hbefore hafter)

theorem HenselTreePreservesFrom.trans
    {lower : Nat} {first second third : Array HenselNode}
    (hfirst : HenselTreePreservesFrom lower first second)
    (hsecond : HenselTreePreservesFrom lower second third) :
    HenselTreePreservesFrom lower first third := by
  refine ⟨le_trans hfirst.1 hsecond.1, ?_⟩
  intro index hlower hindex
  exact (hsecond.2 index hlower (lt_of_lt_of_le hindex hfirst.1)).trans
    (hfirst.2 index hlower hindex)

/-- Constructor-shaped certificate emitted alongside the concrete preorder
builder.  It records raw source pointers; their signed interpretation is
deliberately postponed until the entry-level `h_fits_int32` premise. -/
inductive HenselTreeRawBuildCertificate (lower : Nat) :
    Nat → Nat → Nat → Array HenselNode → Prop
  | node
      (start stop parent : Nat) (nodes : Array HenselNode)
      (value : HenselNode)
      (hlower : lower ≤ parent)
      (hnode : nodes[parent]? = some value)
      (htopology : HenselTreeNodeRawTopologyMatches start stop parent value)
      (hleft : 2 ≤ (start + stop) / 2 - start →
        HenselTreeRawBuildCertificate lower start ((start + stop) / 2)
          (parent + 1) nodes)
      (hright : 2 ≤ stop - (start + stop) / 2 →
        HenselTreeRawBuildCertificate lower ((start + stop) / 2) stop
          (parent + 1 + henselTreeInternalNodeCount start
            ((start + stop) / 2)) nodes) :
      HenselTreeRawBuildCertificate lower start stop parent nodes

/-- Constructor-shaped algebraic companion to the raw topology certificate.
Every node is tied to the interval products computed by the actual builder,
and recursive premises use exactly the two source half-intervals. -/
inductive HenselTreeSemanticBuildCertificate (p : Nat) [Fact (Nat.Prime p)]
    (factors : Array SparsePolyZp) (lower : Nat) :
    (start stop : Nat) → Generated.StrictHensel.HenselLiftTree →
      Array HenselNode → Prop
  | node
      (start stop index : Nat)
      (left right : Option Generated.StrictHensel.HenselLiftTree)
      (nodes : Array HenselNode) (value : HenselNode)
      (hlower : lower ≤ index)
      (hnode : nodes[index]? = some value)
      (hinvariant : HenselTreeNodeGCDInvariant p factors start stop value)
      (hOneHead : HasPhysicalOneHead value.h)
      (hleftCertificate : ∀ child, left = some child →
        HenselTreeSemanticBuildCertificate p factors lower start
          ((start + stop) / 2) child nodes)
      (hrightCertificate : ∀ child, right = some child →
        HenselTreeSemanticBuildCertificate p factors lower
          ((start + stop) / 2) stop child nodes) :
      HenselTreeSemanticBuildCertificate p factors lower start stop
        (.node index left right) nodes

theorem HenselTreeSemanticBuildCertificate.lower_mono
    {p lower lower' start stop : Nat} [Fact (Nat.Prime p)]
    {factors : Array SparsePolyZp}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {nodes : Array HenselNode} (hlower : lower' ≤ lower)
    (hcertificate : HenselTreeSemanticBuildCertificate p factors lower
      start stop tree nodes) :
    HenselTreeSemanticBuildCertificate p factors lower' start stop tree
      nodes := by
  induction hcertificate with
  | node start stop index left right nodes value hindex hnode hinvariant hOneHead
      hleftCertificate hrightCertificate leftIH rightIH =>
      exact .node start stop index left right nodes value (by omega) hnode
        hinvariant hOneHead (fun child hchild => leftIH child hchild)
        (fun child hchild => rightIH child hchild)

theorem henselFactorRangeProduct_singleton
    (p : Nat) (factors : Array SparsePolyZp) (start stop : Nat)
    (hlength : stop - start = 1) :
    henselFactorRangeProduct p factors stop start =
      CLPoly.Math.SparsePolyZp.toPoly p factors[start]! := by
  have hstop : stop = start + 1 := by omega
  subst stop
  rw [henselFactorRangeProduct, if_pos (by omega),
    henselFactorRangeProduct, if_neg (by omega)]
  simp

/-- Every leaf read from the concrete canonical builder array comes from the
corresponding adjusted finite-field input factor, in exact source order. -/
theorem HenselTreeSemanticBuildCertificate.extractedFactors_forall₂
    {p lower start stop root : Nat} [Fact (Nat.Prime p)]
    {factors : Array SparsePolyZp} {nodes : Array HenselNode}
    (hcertificate : HenselTreeSemanticBuildCertificate p factors lower
      start stop (henselTreeBuildTopology start stop root) nodes)
    (hlength : 2 ≤ stop - start) :
    List.Forall₂
      (fun input lifted => toPolyMod p lifted =
        CLPoly.Math.SparsePolyZp.toPoly p input)
      (henselFactorRangeList factors stop start)
      (henselExtractedFactors
        (henselTreeBuildTopology start stop root) nodes) := by
  let mid := (start + stop) / 2
  have hstartMid : start < mid := by
    dsimp [mid]
    exact Generated.StrictHensel.henselTreeMidpoint_gt_start start stop hlength
  have hmidStop : mid < stop := by
    dsimp [mid]
    exact Generated.StrictHensel.henselTreeMidpoint_lt_stop start stop hlength
  have hsplit := henselFactorRangeList_split factors start mid stop
    (by omega) (by omega)
  by_cases hleft : 2 ≤ mid - start
  · by_cases hright : 2 ≤ stop - mid
    · rw [henselTreeBuildTopology] at hcertificate ⊢
      rw [if_pos hleft, if_pos hright] at hcertificate ⊢
      cases hcertificate with
      | node _ _ _ _ _ _ value _ hnode hinvariant _ hleftCert hrightCert =>
          rw [henselExtractedFactors.eq_def]
          rw [hsplit]
          apply List.rel_append
          · exact extractedFactors_forall₂ (hleftCert _ rfl) hleft
          · exact extractedFactors_forall₂ (hrightCert _ rfl) hright
    · have hrightLength : stop - mid = 1 := by omega
      rw [henselTreeBuildTopology] at hcertificate ⊢
      rw [if_pos hleft, if_neg hright] at hcertificate ⊢
      cases hcertificate with
      | node _ _ _ _ _ _ value _ hnode hinvariant _ hleftCert hrightCert =>
          have hindex : root < nodes.size := by
            by_contra hnot
            rw [Array.getElem?_eq_none (by omega)] at hnode
            contradiction
          have hlookup : nodes[root]! = value := by
            have heq : nodes[root] = value :=
              Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
            simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
          rw [henselExtractedFactors.eq_def]
          rw [hsplit, henselFactorRangeList_singleton factors mid stop
            hrightLength]
          apply List.rel_append
          · exact extractedFactors_forall₂ (hleftCert _ rfl) hleft
          · exact .cons (by
              rw [hlookup]
              exact hinvariant.2.1.trans
                (henselFactorRangeProduct_singleton p factors mid stop
                  hrightLength)) .nil
  · have hleftLength : mid - start = 1 := by omega
    by_cases hright : 2 ≤ stop - mid
    · rw [henselTreeBuildTopology] at hcertificate ⊢
      rw [if_neg hleft, if_pos hright] at hcertificate ⊢
      cases hcertificate with
      | node _ _ _ _ _ _ value _ hnode hinvariant _ hleftCert hrightCert =>
          have hindex : root < nodes.size := by
            by_contra hnot
            rw [Array.getElem?_eq_none (by omega)] at hnode
            contradiction
          have hlookup : nodes[root]! = value := by
            have heq : nodes[root] = value :=
              Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
            simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
          rw [henselExtractedFactors.eq_def]
          rw [hsplit, henselFactorRangeList_singleton factors start mid
            hleftLength]
          apply List.rel_append
          · exact .cons (by
              rw [hlookup]
              exact hinvariant.1.trans
                (henselFactorRangeProduct_singleton p factors start mid
                  hleftLength)) .nil
          · exact extractedFactors_forall₂ (hrightCert _ rfl) hright
    · have hrightLength : stop - mid = 1 := by omega
      rw [henselTreeBuildTopology] at hcertificate ⊢
      rw [if_neg hleft, if_neg hright] at hcertificate ⊢
      cases hcertificate with
      | node _ _ _ _ _ _ value _ hnode hinvariant _ hleftCert hrightCert =>
          have hindex : root < nodes.size := by
            by_contra hnot
            rw [Array.getElem?_eq_none (by omega)] at hnode
            contradiction
          have hlookup : nodes[root]! = value := by
            have heq : nodes[root] = value :=
              Option.some.inj ((Array.getElem?_eq_getElem hindex).symm.trans hnode)
            simpa [getElem!_def, Array.getElem?_eq_getElem hindex] using heq
          rw [henselExtractedFactors.eq_def]
          rw [hsplit,
            henselFactorRangeList_singleton factors start mid hleftLength,
            henselFactorRangeList_singleton factors mid stop hrightLength]
          exact .cons (by
              rw [hlookup]
              exact hinvariant.1.trans
                (henselFactorRangeProduct_singleton p factors start mid
                  hleftLength))
            (.cons (by
              rw [hlookup]
              exact hinvariant.2.1.trans
                (henselFactorRangeProduct_singleton p factors mid stop
                  hrightLength)) .nil)
termination_by stop - start
decreasing_by all_goals omega

/-- Pointwise leaf origins preserve the ordered product.  This merely folds
the exact `Forall₂` correspondence; it does not factor or reconstruct a
polynomial. -/
private theorem henselLeafOrigins_product
    {p : Nat} {inputs : List SparsePolyZp} {lifted : List SparsePolyZZ}
    (horigins : List.Forall₂
      (fun input output => toPolyMod p output =
        CLPoly.Math.SparsePolyZp.toPoly p input) inputs lifted) :
    (lifted.map (toPolyMod p)).prod =
      (inputs.map (CLPoly.Math.SparsePolyZp.toPoly p)).prod := by
  induction horigins with
  | nil => simp
  | cons hhead htail ih => simp [hhead, ih]

/-- The leaves stored by the actual generated builder multiply to the product
of the exact source factor interval, in source order. -/
theorem HenselTreeSemanticBuildCertificate.extractedFactors_product
    {p lower start stop root : Nat} [Fact (Nat.Prime p)]
    {factors : Array SparsePolyZp} {nodes : Array HenselNode}
    (hcertificate : HenselTreeSemanticBuildCertificate p factors lower
      start stop (henselTreeBuildTopology start stop root) nodes)
    (hlength : 2 ≤ stop - start) :
    ((henselExtractedFactors (henselTreeBuildTopology start stop root)
      nodes).map (toPolyMod p)).prod =
      ((henselFactorRangeList factors stop start).map
        (CLPoly.Math.SparsePolyZp.toPoly p)).prod := by
  exact henselLeafOrigins_product
    (hcertificate.extractedFactors_forall₂ hlength)

/-- Combining the concrete builder certificate with the concrete
well-founded quadratic loop transports a supplied source-factor product
equation to the exact modulus and node array returned by that loop.  The
source equation is deliberately explicit: it is discharged by the preceding
finite-field factorization refinement, not postulated by Hensel lifting. -/
theorem HenselTreeSemanticBuildCertificate.liftLoop_extractedFactors_product
    {termination : Generated.StrictHensel.DivmodTermination}
    {p lower start stop root target outputM : Nat} [Fact (Nat.Prime p)]
    {factors : Array SparsePolyZp} {f : SparsePolyZZ}
    {nodes outputNodes : Array HenselNode}
    (hcertificate : HenselTreeSemanticBuildCertificate p factors lower
      start stop (henselTreeBuildTopology start stop root) nodes)
    (hloop : HenselLiftLoopCorrect termination
      (henselTreeBuildTopology start stop root) f target p nodes outputNodes
      outputM)
    (hlength : 2 ≤ stop - start)
    (hsource : ((henselFactorRangeList factors stop start).map
      (CLPoly.Math.SparsePolyZp.toPoly p)).prod = toPolyMod p f) :
    ((henselExtractedFactors (henselTreeBuildTopology start stop root)
      outputNodes).map (toPolyMod outputM)).prod =
        toPolyMod outputM f := by
  apply hloop.extractedFactors_product
  · exact (henselTreeBuildTopology_indices_nodup_bounded start stop root
      hlength).1
  · exact (hcertificate.extractedFactors_product hlength).trans hsource

/-- Reading a certified builder topology before and after the actual Hensel
lift produces pointwise congruent leaf lists modulo every divisor tracked by
`HenselArrayReduces`. -/
theorem HenselTreeSemanticBuildCertificate.extractedFactors_forall₂_of_reduces
    {p lower start stop : Nat} [Fact (Nat.Prime p)]
    {factors : Array SparsePolyZp}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (hcertificate : HenselTreeSemanticBuildCertificate p factors lower
      start stop tree before)
    (hreduce : HenselArrayReduces p before after) :
    List.Forall₂ (fun old new => toPolyMod p new = toPolyMod p old)
      (henselExtractedFactors tree before)
      (henselExtractedFactors tree after) := by
  induction hcertificate with
  | node start stop index left right nodes value hlower hnode hinvariant hOneHead
      hleftCertificate hrightCertificate leftIH rightIH =>
      have hindex : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      have hafterIndex : index < after.size := by
        rw [← hreduce.1]
        exact hindex
      have hnodeReduce := hreduce.2 index hindex hafterIndex
      have hnodeReduceBang :
          HenselNodeReduces p nodes[index]! after[index]! := by
        simpa [getElem!_def, Array.getElem?_eq_getElem hindex,
          Array.getElem?_eq_getElem hafterIndex] using hnodeReduce
      rw [henselExtractedFactors.eq_def, henselExtractedFactors.eq_def]
      apply List.rel_append
      · cases left with
        | none => exact .cons hnodeReduceBang.1 .nil
        | some child => exact leftIH child rfl hreduce
      · cases right with
        | none => exact .cons hnodeReduceBang.2.1 .nil
        | some child => exact rightIH child rfl hreduce

theorem HenselTreeSemanticBuildCertificate.of_preservesFrom
    {p lower start stop : Nat} [Fact (Nat.Prime p)]
    {factors : Array SparsePolyZp}
    {tree : Generated.StrictHensel.HenselLiftTree}
    {before after : Array HenselNode}
    (hpreserves : HenselTreePreservesFrom lower before after)
    (hcertificate : HenselTreeSemanticBuildCertificate p factors lower
      start stop tree before) :
    HenselTreeSemanticBuildCertificate p factors lower start stop tree
      after := by
  induction hcertificate with
  | node start stop index left right nodes value hindex hnode hinvariant hOneHead
      hleftCertificate hrightCertificate leftIH rightIH =>
      have hbound : index < nodes.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hnode
        contradiction
      exact .node start stop index left right after value hindex
        ((hpreserves.2 index hindex hbound).trans hnode) hinvariant hOneHead
        (fun child hchild => leftIH child hchild hpreserves)
        (fun child hchild => rightIH child hchild hpreserves)

theorem HenselTreeRawBuildCertificate.lower_mono
    {lower lower' start stop parent : Nat} {nodes : Array HenselNode}
    (hlower : lower' ≤ lower)
    (hcertificate : HenselTreeRawBuildCertificate lower start stop parent nodes) :
    HenselTreeRawBuildCertificate lower' start stop parent nodes := by
  induction hcertificate with
  | node start stop parent nodes value hparent hnode htopology hleft hright
      leftIH rightIH =>
      exact .node start stop parent nodes value (by omega) hnode htopology
        (fun h => leftIH h) (fun h => rightIH h)

theorem HenselTreeRawBuildCertificate.of_preservesFrom
    {lower start stop parent : Nat} {before after : Array HenselNode}
    (hpreserves : HenselTreePreservesFrom lower before after)
    (hcertificate :
      HenselTreeRawBuildCertificate lower start stop parent before) :
    HenselTreeRawBuildCertificate lower start stop parent after := by
  induction hcertificate with
  | node start stop parent nodes value hlower hnode htopology hleft hright
      leftIH rightIH =>
      have hparent : parent < nodes.size := by
        by_contra hnot
        have : nodes[parent]? = none := Array.getElem?_eq_none (by omega)
        rw [this] at hnode
        contradiction
      exact .node start stop parent after value hlower
        ((hpreserves.2 parent hlower hparent).trans hnode) htopology
        (fun h => leftIH h hpreserves) (fun h => rightIH h hpreserves)

set_option maxHeartbeats 0 in
/-- Interpreting every raw pointer is valid once the whole preorder allocation
fits in signed `Int32`; the result is exactly the certificate consumed by the
generated extraction traversal. -/
theorem HenselTreeRawBuildCertificate.toExtractCertificate
    {lower start stop parent : Nat} {nodes : Array HenselNode}
    (hlength : 2 ≤ stop - start)
    (hbound : parent + henselTreeInternalNodeCount start stop < 2 ^ 31)
    (hcertificate :
      HenselTreeRawBuildCertificate lower start stop parent nodes) :
    HenselExtractCertificate lower
      (henselTreeBuildTopology start stop parent) nodes := by
  cases hcertificate with
  | node _ _ _ _ value hlower hnode htopology hleft hright =>
      have hchildren := htopology.liftChildMatches hlength hbound
      let mid := (start + stop) / 2
      by_cases hleftBranch : 2 ≤ mid - start
      · have hleftCount := henselTreeInternalNodeCount_eq start mid hleftBranch
        by_cases hrightBranch : 2 ≤ stop - mid
        · have hrightCount := henselTreeInternalNodeCount_eq mid stop hrightBranch
          dsimp [mid] at hleftBranch hrightBranch
          rw [henselTreeBuildTopology] at hchildren ⊢
          dsimp only at hchildren ⊢
          rw [if_pos hleftBranch, if_pos hrightBranch] at hchildren ⊢
          refine .node parent _ _ nodes value hlower hnode hchildren.1
            hchildren.2 ?_ ?_
          · intro child hchild
            injection hchild with hchildEq
            subst child
            exact (hleft hleftBranch).toExtractCertificate hleftBranch (by
              have htotal := henselTreeInternalNodeCount_eq start stop hlength
              dsimp [mid] at htotal hbound ⊢
              omega)
          · intro child hchild
            injection hchild with hchildEq
            subst child
            exact (hright hrightBranch).toExtractCertificate hrightBranch (by
              have htotal := henselTreeInternalNodeCount_eq start stop hlength
              dsimp [mid] at htotal hbound ⊢
              omega)
        · dsimp [mid] at hleftBranch hrightBranch
          rw [henselTreeBuildTopology] at hchildren ⊢
          dsimp only at hchildren ⊢
          rw [if_pos hleftBranch, if_neg hrightBranch] at hchildren ⊢
          refine .node parent _ _ nodes value hlower hnode hchildren.1
            hchildren.2 ?_ ?_
          · intro child hchild
            injection hchild with hchildEq
            subst child
            exact (hleft hleftBranch).toExtractCertificate hleftBranch (by
              have htotal := henselTreeInternalNodeCount_eq start stop hlength
              dsimp [mid] at htotal hbound ⊢
              omega)
          · intro child hchild
            contradiction
      · by_cases hrightBranch : 2 ≤ stop - mid
        · dsimp [mid] at hleftBranch hrightBranch
          rw [henselTreeBuildTopology] at hchildren ⊢
          dsimp only at hchildren ⊢
          rw [if_neg hleftBranch, if_pos hrightBranch] at hchildren ⊢
          refine .node parent _ _ nodes value hlower hnode hchildren.1
            hchildren.2 ?_ ?_
          · intro child hchild
            contradiction
          · intro child hchild
            injection hchild with hchildEq
            subst child
            exact (hright hrightBranch).toExtractCertificate hrightBranch (by
              have htotal := henselTreeInternalNodeCount_eq start stop hlength
              have hleftZero : henselTreeInternalNodeCount start mid = 0 := by
                rw [henselTreeInternalNodeCount, dif_neg hleftBranch]
              dsimp [mid] at htotal hleftZero hbound ⊢
              omega)
        · dsimp [mid] at hleftBranch hrightBranch
          rw [henselTreeBuildTopology] at hchildren ⊢
          dsimp only at hchildren ⊢
          rw [if_neg hleftBranch, if_neg hrightBranch] at hchildren ⊢
          refine .node parent _ _ nodes value hlower hnode hchildren.1
            hchildren.2 ?_ ?_
          · intro child hchild
            contradiction
          · intro child hchild
            contradiction

theorem henselTreeSetNodeRawIR_frame_current
    (nodes : Array HenselNode) (index : Nat)
    (update : HenselNode → HenselNode) (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeSetNodeRawIR nodes index update =
        .ok output ∧
      HenselTreeArrayFrame nodes output index := by
  rcases henselTreeSetNodeRawIR_frame nodes index update hindex with
    ⟨output, hrun, hsize, houtput, hselected, hother⟩
  refine ⟨output, hrun, ⟨by omega, ?_⟩⟩
  intro other hne hbefore hafter
  exact hother other hbefore hne hafter

theorem henselTreeSetNodeRawIR_preserves_algebra
    (nodes : Array HenselNode) (index : Nat)
    (update : HenselNode → HenselNode) (hindex : index < nodes.size)
    (hupdate : HenselTreeNodeAlgebraEq (update nodes[index]) nodes[index]) :
    ∃ output,
      Generated.StrictHensel.henselTreeSetNodeRawIR nodes index update =
        .ok output ∧
      output.size = nodes.size ∧
      HenselTreeArrayFrame nodes output index ∧
      ∃ houtput : index < output.size,
      HenselTreeNodeAlgebraEq (getElem output index houtput) nodes[index] ∧
      getElem output index houtput = update nodes[index] := by
  rcases henselTreeSetNodeRawIR_frame nodes index update hindex with
    ⟨output, hrun, hsize, houtputIndex, hselected, hother⟩
  refine ⟨output, hrun, hsize, ⟨by omega, ?_⟩, houtputIndex, ?_, ?_⟩
  · intro other hne hbefore hafter
    exact hother other hbefore hne hafter
  · simpa only [hselected] using hupdate
  · exact hselected

theorem henselTreeStoreNodeRawIR_frame_current
    (nodes : Array HenselNode) (index : Nat)
    (g h s t : SparsePolyZp) (start stop : Nat)
    (hindex : index < nodes.size) :
    ∃ output,
      Generated.StrictHensel.henselTreeStoreNodeRawIR nodes index g h s t
          start stop = .ok output ∧
      HenselTreeArrayFrame nodes output index := by
  rcases henselTreeStoreNodeRawIR_frame nodes index g h s t start stop
      hindex with ⟨output, hrun, hsize, hother⟩
  refine ⟨output, hrun, ⟨by omega, ?_⟩⟩
  intro other hne hbefore hafter
  exact hother other hbefore hne hafter

theorem henselTreeArrayFrame_push (nodes : Array HenselNode)
    (current : Nat) :
    HenselTreeArrayFrame nodes (nodes.push default) current := by
  refine ⟨by simp, ?_⟩
  intro index hne hbefore hafter
  simp only [Array.getElem_push, hbefore, ↓reduceDIte]

/-- Execution-only safety for the concrete well-founded tree builder.  All
products and EEA calls are the actual strict operations, and recursive calls
can only append nodes. -/
theorem strictHenselTreeBuildRecursiveRawIR_succeeds
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic) :
    ∀ nodes start stop parent,
      2 ≤ stop - start → stop ≤ factors.size → parent < nodes.size →
      nodes.size = parent + 1 → HenselArrayCanonical nodes →
      ∃ output,
        Generated.StrictHensel.__hensel_tree_build_recursive_raw_ir
            (strictHenselTreeBuildRawOps this mulProvider) factors nodes start
            stop parent = .ok output ∧
        HenselTreeArrayFrame nodes output parent ∧
        ∃ hparentOutput : parent < output.size,
        output.size = nodes.size +
          henselTreeInternalNodeCount start ((start + stop) / 2) +
          henselTreeInternalNodeCount ((start + stop) / 2) stop ∧
        HenselTreeNodeRawTopologyMatches start stop parent
          (getElem output parent hparentOutput) ∧
        HenselTreeNodeGCDInvariant this._p.toNat factors start stop
          (getElem output parent hparentOutput) ∧
        HenselTreeRawBuildCertificate parent start stop parent output ∧
        HenselTreeSemanticBuildCertificate this._p.toNat factors parent
          start stop (henselTreeBuildTopology start stop parent) output ∧
        HenselArrayCanonical output := by
  intro nodes start stop parent hlength hstop hparent hparentLast hnodesCanonical
  rw [Generated.StrictHensel.__hensel_tree_build_recursive_raw_ir]
  rw [dif_pos hlength]
  let mid := (start + stop) / 2
  dsimp only
  have hstartMid : start < mid :=
    Generated.StrictHensel.henselTreeMidpoint_gt_start start stop hlength
  have hmidStop : mid < stop :=
    Generated.StrictHensel.henselTreeMidpoint_lt_stop start stop hlength
  rcases strictHenselTreeProductRangeRawIR_refines this hcfg mulProvider
      factors hfactors hfactorsNonempty start mid hstartMid (by omega) with
    ⟨g, hgRun, hgCanonical, hgNonempty, hgPoly⟩
  rcases strictHenselTreeProductRangeRawIR_refines this hcfg mulProvider
      factors hfactors hfactorsNonempty mid stop hmidStop hstop with
    ⟨h, hhRun, hhCanonical, hhNonempty, hhPoly⟩
  have hhMonic :
      (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat h).Monic := by
    rw [hhPoly]
    exact henselFactorRangeProduct_monic this._p.toNat factors mid stop hstop
      (by
        intro index hindex hmidIndex hindexStop
        exact hfactorsMonicAfterZero index hindex (by omega))
  rw [hgRun, hhRun]
  rcases strictHenselEEAEntryIR_refines_gcd this hcfg h2p hp2 g h
      hgCanonical hhCanonical hgNonempty with
    ⟨gcd, s, t, heeaRun, hbezout, hmonic, hdvdG, hdvdH, hgcd⟩
  simp only [strictHenselTreeBuildRawOps, heeaRun, bind, Except.bind]
  rcases henselTreeStoreNodeRawIR_refines_gcd this._p.toNat factors nodes
      parent start stop g h gcd s t hparent (by simpa [mid] using hgPoly)
      (by simpa [mid] using hhPoly) hbezout hmonic hgcd with
    ⟨stored, hstoredRun, hparentStored, hstoredSize, hstoredInvariant⟩
  rcases henselTreeStoreNodeRawIR_frame_current nodes parent g h s t start stop
      hparent with ⟨stored', hstoredRun', hstoredFrame⟩
  rw [hstoredRun] at hstoredRun'
  injection hstoredRun' with hstoredEq
  subst stored'
  rcases henselTreeStoreNodeRawIR_canonical_oneHead this._p.toNat nodes parent
      g h s t start stop hparent hgCanonical hhCanonical hhNonempty hhMonic with
    ⟨stored', hstoredCanonicalRun, hstoredCanonicalSize,
      hparentStored', hstoredCanonical, hstoredOneHead⟩
  rw [hstoredRun] at hstoredCanonicalRun
  injection hstoredCanonicalRun with hstoredCanonicalEq
  subst stored'
  have hstoredArrayCanonical : HenselArrayCanonical stored :=
    hstoredFrame.canonical_of_size_eq hstoredCanonicalSize
      hnodesCanonical hparentStored' hstoredCanonical
  simp only [hstoredRun, bind, Except.bind]
  by_cases hleft : 2 ≤ (start + stop) / 2 - start
  · rw [dif_pos hleft]
    let child := stored.size
    let pushed := stored.push default
    have hframePushed : HenselTreeArrayFrame nodes pushed parent :=
      hstoredFrame.trans (henselTreeArrayFrame_push stored parent)
    have hparentPushed : parent < pushed.size := by
      simp [pushed]
      omega
    rcases henselTreeSetNodeRawIR_preserves_algebra pushed parent
        (fun node => { node with left := child.toUInt32.toInt32 })
        hparentPushed (by simp [HenselTreeNodeAlgebraEq]) with
      ⟨leftReady, hleftReadyRun, hleftReadySizeFromRun, hleftReadyFrame,
        hleftReadyParent,
        hleftReadyAlgebra, hleftReadySelected⟩
    have hframeLeftReady : HenselTreeArrayFrame nodes leftReady parent :=
      hframePushed.trans hleftReadyFrame
    have hleftReadySize : leftReady.size = pushed.size :=
      hleftReadySizeFromRun
    have hpushedParent : pushed[parent] = stored[parent] := by
      simp only [pushed, Array.getElem_push, hparentStored, ↓reduceDIte]
    have hleftReadyStored :
        HenselTreeNodeAlgebraEq leftReady[parent] stored[parent] :=
      hleftReadyAlgebra.trans (by simpa [hpushedParent] using
        (HenselTreeNodeAlgebraEq.refl stored[parent]))
    have hpushedCanonical : HenselArrayCanonical pushed :=
      hstoredArrayCanonical.push_default
    have hleftReadyCanonical : HenselArrayCanonical leftReady :=
      hleftReadyFrame.canonical_of_size_eq hleftReadySize hpushedCanonical
        hleftReadyParent (hstoredCanonical.of_algebraEq hleftReadyStored)
    dsimp [pushed, child] at hleftReadyRun
    rw [hleftReadyRun]
    have hchildLeftReady : child < leftReady.size := by
      simp [child, pushed] at hleftReadySize ⊢
      omega
    rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
        mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
        leftReady start mid child
        hleft (by omega) hchildLeftReady (by
          simp [child, pushed] at hleftReadySize ⊢
          omega) hleftReadyCanonical with
      ⟨afterLeft, hafterLeftRun, hafterLeftFrame, hparentAfterLeft,
        hafterLeftSizeExact, hafterLeftTopology, hafterLeftInvariant,
        hafterLeftCertificate, hafterLeftSemantic, hafterLeftCanonical⟩
    have hframeAfterLeft : HenselTreeArrayFrame nodes afterLeft parent :=
      hframeLeftReady.trans_fresh hafterLeftFrame (by
        dsimp [child]
        omega)
    have hleftCountEq := henselTreeInternalNodeCount_eq start mid hleft
    have hparentNeChild : parent ≠ child := by
      dsimp [child]
      omega
    have hafterLeftParent : afterLeft[parent] = leftReady[parent] :=
      hafterLeftFrame.2 parent hparentNeChild (by omega) (by omega)
    have hafterLeftStored :
        HenselTreeNodeAlgebraEq afterLeft[parent] stored[parent] :=
      (by simpa [hafterLeftParent] using hleftReadyStored)
    simp only [bind, Except.bind]
    have hafterLeftRunExpanded :
      Generated.StrictHensel.__hensel_tree_build_recursive_raw_ir
          { mul := fun left right => StrictDDF.strictMulIR this left right mulProvider
            eea := fun left right =>
              Generated.StrictHensel.__polynomial_GCD_eea_raw_ir
                (strictHenselEEARawOps this)
                (strictHenselEEATermination this)
                (Generated.StrictHensel.henselEEAInitialState this._p left right) }
          factors leftReady start ((start + stop) / 2) stored.size =
            .ok afterLeft := by
      simpa [strictHenselTreeBuildRawOps, mid, child] using hafterLeftRun
    rw [hafterLeftRunExpanded]
    simp only [bind, Except.bind]
    by_cases hright : 2 ≤ stop - (start + stop) / 2
    · rw [dif_pos hright]
      let rightChild := afterLeft.size
      let rightPushed := afterLeft.push default
      have hframeRightPushed : HenselTreeArrayFrame nodes rightPushed parent :=
        hframeAfterLeft.trans (henselTreeArrayFrame_push afterLeft parent)
      have hparentRightPushed : parent < rightPushed.size := by
        simp [rightPushed]
        omega
      rcases henselTreeSetNodeRawIR_preserves_algebra rightPushed parent
          (fun node => { node with right := rightChild.toUInt32.toInt32 })
          hparentRightPushed (by simp [HenselTreeNodeAlgebraEq]) with
        ⟨rightReady, hrightReadyRun, hrightReadySizeFromRun,
          hrightReadyFrame, hrightReadyParent,
          hrightReadyAlgebra, hrightReadySelected⟩
      have hframeRightReady : HenselTreeArrayFrame nodes rightReady parent :=
        hframeRightPushed.trans hrightReadyFrame
      have hrightReadySize : rightReady.size = rightPushed.size :=
        hrightReadySizeFromRun
      have hrightPushedParent : rightPushed[parent] = afterLeft[parent] := by
        dsimp only [rightPushed]
        exact @Array.getElem_push_lt HenselNode afterLeft default parent
          (by omega)
      have hrightReadyStored :
          HenselTreeNodeAlgebraEq rightReady[parent] stored[parent] :=
        hrightReadyAlgebra.trans
          ((HenselTreeNodeAlgebraEq.of_eq hrightPushedParent).trans
            hafterLeftStored)
      have hrightPushedCanonical : HenselArrayCanonical rightPushed :=
        hafterLeftCanonical.push_default
      have hrightReadyCanonical : HenselArrayCanonical rightReady :=
        hrightReadyFrame.canonical_of_size_eq hrightReadySize
          hrightPushedCanonical hrightReadyParent
          (hstoredCanonical.of_algebraEq hrightReadyStored)
      dsimp [rightPushed, rightChild] at hrightReadyRun
      rw [hrightReadyRun]
      have hrightChildReady : rightChild < rightReady.size := by
        simp [rightChild, rightPushed] at hrightReadySize ⊢
        omega
      rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
          mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
          rightReady mid stop
          rightChild hright hstop hrightChildReady (by
            simp [rightChild, rightPushed] at hrightReadySize ⊢
            omega) hrightReadyCanonical with
        ⟨output, houtputRun, houtputFrame, hparentOutputBound,
          houtputSizeExact, houtputTopology, houtputInvariant,
          houtputCertificate, houtputSemantic, houtputCanonical⟩
      have hrightCountEq := henselTreeInternalNodeCount_eq mid stop hright
      have hparentNeRightChild : parent ≠ rightChild := by
        dsimp [rightChild]
        omega
      have houtputParent : output[parent] = rightReady[parent] :=
        houtputFrame.2 parent hparentNeRightChild (by omega) (by omega)
      exact ⟨output, houtputRun,
        hframeRightReady.trans_fresh houtputFrame (by
          dsimp [rightChild]
          omega),
        by omega,
        by
          dsimp [mid] at hleftCountEq hrightCountEq hafterLeftSizeExact houtputSizeExact ⊢
          simp [rightPushed] at hrightReadySize
          simp [pushed] at hleftReadySize
          omega,
        by
          have hchildEq : child = parent + 1 := by
            dsimp [child]
            omega
          have hrightChildEq :
              rightChild = parent + 1 +
                henselTreeInternalNodeCount start ((start + stop) / 2) := by
            dsimp [rightChild, mid] at hleftCountEq hafterLeftSizeExact ⊢
            simp [pushed] at hleftReadySize
            omega
          simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
            houtputParent, hrightReadySelected, hrightPushedParent,
            hafterLeftParent, hleftReadySelected, hchildEq, hrightChildEq],
        hstoredInvariant.of_algebraEq (by
          simpa [houtputParent] using hrightReadyStored),
        by
          have hafterLeftOutputFrame :
              HenselTreeArrayFrame afterLeft output parent :=
            ((henselTreeArrayFrame_push afterLeft parent).trans
              hrightReadyFrame).trans_fresh houtputFrame (by
                dsimp [rightChild]
                omega)
          have hleftFinal := hafterLeftCertificate.of_preservesFrom
            (hafterLeftOutputFrame.preservesFrom (by
              dsimp [child]
              omega))
          have hchildEq : child = parent + 1 := by
            dsimp [child]
            omega
          have hrightChildEq : rightChild = parent + 1 +
              henselTreeInternalNodeCount start ((start + stop) / 2) := by
            dsimp [rightChild, mid] at hleftCountEq hafterLeftSizeExact ⊢
            simp [pushed] at hleftReadySize
            omega
          have hleftFinal' := hleftFinal.lower_mono (by omega : parent ≤ child)
          have hrightFinal := houtputCertificate.lower_mono
            (by omega : parent ≤ rightChild)
          exact .node start stop parent output output[parent] (by omega)
            (@Array.getElem?_eq_getElem _ output parent (by omega))
            (by
              simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
                houtputParent, hrightReadySelected, hrightPushedParent,
                hafterLeftParent, hleftReadySelected, hchildEq,
                hrightChildEq])
            (fun _ => by simpa [hchildEq] using hleftFinal')
            (fun _ => by simpa [hrightChildEq] using hrightFinal),
        by
          have hafterLeftOutputFrame :
              HenselTreeArrayFrame afterLeft output parent :=
            ((henselTreeArrayFrame_push afterLeft parent).trans
              hrightReadyFrame).trans_fresh houtputFrame (by
                dsimp [rightChild]
                omega)
          have hleftFinal := hafterLeftSemantic.of_preservesFrom
            (hafterLeftOutputFrame.preservesFrom (by
              dsimp [child]
              omega))
          have hchildEq : child = parent + 1 := by
            dsimp [child]
            omega
          have hrightChildEq : rightChild = parent + 1 +
              henselTreeInternalNodeCount start ((start + stop) / 2) := by
            dsimp [rightChild, mid] at hleftCountEq hafterLeftSizeExact ⊢
            simp [pushed] at hleftReadySize
            omega
          have hleftFinal' := hleftFinal.lower_mono
            (by omega : parent ≤ child)
          have hrightFinal := houtputSemantic.lower_mono
            (by omega : parent ≤ rightChild)
          rw [henselTreeBuildTopology]
          rw [if_pos hleft, if_pos hright]
          exact .node start stop parent _ _ output output[parent] (by omega)
            (@Array.getElem?_eq_getElem _ output parent (by omega))
            (hstoredInvariant.of_algebraEq (by
              simpa [houtputParent] using hrightReadyStored))
            (hstoredOneHead.of_algebraEq (by
              simpa [houtputParent] using hrightReadyStored))
            (fun childTree hchildTree => by
              injection hchildTree with htree
              subst childTree
              simpa [hchildEq] using hleftFinal')
            (fun childTree hchildTree => by
              injection hchildTree with htree
              subst childTree
              simpa [hrightChildEq] using hrightFinal),
        houtputCanonical⟩
    · rw [dif_neg hright]
      have hparentAfterLeft : parent < afterLeft.size := by omega
      rcases henselTreeSetNodeRawIR_preserves_algebra afterLeft parent
          (fun node => { node with right := -1 }) hparentAfterLeft
          (by simp [HenselTreeNodeAlgebraEq]) with
        ⟨output, houtputRun, houtputSizeFromRun, houtputFrame,
          houtputParentBound,
          houtputAlgebra, houtputSelected⟩
      have houtputCanonical : HenselArrayCanonical output :=
        houtputFrame.canonical_of_size_eq houtputSizeFromRun
          hafterLeftCanonical houtputParentBound
          (hstoredCanonical.of_algebraEq
            (houtputAlgebra.trans hafterLeftStored))
      exact ⟨output, houtputRun, hframeAfterLeft.trans houtputFrame,
        by omega,
        by
          dsimp [mid] at hleftCountEq hafterLeftSizeExact ⊢
          have hrightCount :
              henselTreeInternalNodeCount ((start + stop) / 2) stop = 0 := by
            rw [henselTreeInternalNodeCount, dif_neg hright]
          simp [pushed] at hleftReadySize
          omega,
        by
          simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
            houtputSelected, hafterLeftParent, hleftReadySelected]
          congr
          dsimp [child, mid] at *
          omega,
        hstoredInvariant.of_algebraEq
          (houtputAlgebra.trans hafterLeftStored),
        by
          have hleftFinal := hafterLeftCertificate.of_preservesFrom
            (houtputFrame.preservesFrom (by
              dsimp [child]
              omega))
          have hchildEq : child = parent + 1 := by
            dsimp [child]
            omega
          have hleftFinal' := hleftFinal.lower_mono (by omega : parent ≤ child)
          exact .node start stop parent output output[parent] (by omega)
            (Array.getElem?_eq_getElem houtputParentBound)
            (by
              simpa using (show HenselTreeNodeRawTopologyMatches start stop parent
                output[parent] from by
                  simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
                    houtputSelected, hafterLeftParent, hleftReadySelected]
                  congr
                  dsimp [child, mid] at *
                  omega))
            (fun _ => by simpa [hchildEq] using hleftFinal')
            (fun h => False.elim (hright h)),
        by
          have hleftFinal := hafterLeftSemantic.of_preservesFrom
            (houtputFrame.preservesFrom (by
              dsimp [child]
              omega))
          have hchildEq : child = parent + 1 := by
            dsimp [child]
            omega
          have hleftFinal' := hleftFinal.lower_mono
            (by omega : parent ≤ child)
          rw [henselTreeBuildTopology]
          rw [if_pos hleft, if_neg hright]
          exact .node start stop parent _ _ output output[parent] (by omega)
            (Array.getElem?_eq_getElem houtputParentBound)
            (hstoredInvariant.of_algebraEq
              (houtputAlgebra.trans hafterLeftStored))
            (hstoredOneHead.of_algebraEq
              (houtputAlgebra.trans hafterLeftStored))
            (fun childTree hchildTree => by
              injection hchildTree with htree
              subst childTree
              simpa [hchildEq] using hleftFinal')
            (fun childTree hchildTree => by contradiction),
        houtputCanonical⟩
  · rw [dif_neg hleft]
    have hparentStored : parent < stored.size := by omega
    rcases henselTreeSetNodeRawIR_preserves_algebra stored parent
        (fun node => { node with left := -1 }) hparentStored
        (by simp [HenselTreeNodeAlgebraEq]) with
      ⟨afterLeft, hafterLeftRun, hafterLeftSizeFromRun, hafterLeftFrame,
        hparentAfterLeft,
        hafterLeftAlgebra, hafterLeftSelected⟩
    have hframeAfterLeft : HenselTreeArrayFrame nodes afterLeft parent :=
      hstoredFrame.trans hafterLeftFrame
    have hafterLeftStored :
        HenselTreeNodeAlgebraEq afterLeft[parent] stored[parent] :=
      hafterLeftAlgebra
    have hafterLeftSize : afterLeft.size = stored.size :=
      hafterLeftSizeFromRun
    have hafterLeftCanonical : HenselArrayCanonical afterLeft :=
      hafterLeftFrame.canonical_of_size_eq hafterLeftSize
        hstoredArrayCanonical hparentAfterLeft
        (hstoredCanonical.of_algebraEq hafterLeftStored)
    rw [hafterLeftRun]
    simp only [bind, Except.bind]
    by_cases hright : 2 ≤ stop - (start + stop) / 2
    · rw [dif_pos hright]
      let child := afterLeft.size
      let pushed := afterLeft.push default
      have hframePushed : HenselTreeArrayFrame nodes pushed parent :=
        hframeAfterLeft.trans (henselTreeArrayFrame_push afterLeft parent)
      have hparentPushed : parent < pushed.size := by
        simp [pushed]
        omega
      rcases henselTreeSetNodeRawIR_preserves_algebra pushed parent
          (fun node => { node with right := child.toUInt32.toInt32 })
        hparentPushed (by simp [HenselTreeNodeAlgebraEq]) with
        ⟨rightReady, hrightReadyRun, hrightReadySizeFromRun,
          hrightReadyFrame,
          hrightReadyParent, hrightReadyAlgebra, hrightReadySelected⟩
      have hframeRightReady : HenselTreeArrayFrame nodes rightReady parent :=
        hframePushed.trans hrightReadyFrame
      have hrightReadySizeEq : rightReady.size = pushed.size :=
        hrightReadySizeFromRun
      have hpushedParent : pushed[parent] = afterLeft[parent] := by
        dsimp only [pushed]
        exact @Array.getElem_push_lt HenselNode afterLeft default parent
          hparentAfterLeft
      have hrightReadyStored :
          HenselTreeNodeAlgebraEq rightReady[parent] stored[parent] :=
        hrightReadyAlgebra.trans
          ((HenselTreeNodeAlgebraEq.of_eq hpushedParent).trans
            hafterLeftStored)
      have hpushedCanonical : HenselArrayCanonical pushed :=
        hafterLeftCanonical.push_default
      have hrightReadyCanonical : HenselArrayCanonical rightReady :=
        hrightReadyFrame.canonical_of_size_eq hrightReadySizeEq
          hpushedCanonical hrightReadyParent
          (hstoredCanonical.of_algebraEq hrightReadyStored)
      dsimp [pushed, child] at hrightReadyRun
      rw [hrightReadyRun]
      have hchildReady : child < rightReady.size := by
        simp [child, pushed] at hrightReadySizeEq ⊢
        omega
      rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
          mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
          rightReady mid stop
          child hright hstop hchildReady (by
            simp [child, pushed] at hrightReadySizeEq ⊢
            omega) hrightReadyCanonical with
        ⟨output, houtputRun, houtputFrame, hparentOutputBound,
          houtputSizeExact, houtputTopology, houtputInvariant,
          houtputCertificate, houtputSemantic, houtputCanonical⟩
      have hrightCountEq := henselTreeInternalNodeCount_eq mid stop hright
      have hparentNeChild : parent ≠ child := by
        dsimp [child]
        omega
      have houtputParent : output[parent] = rightReady[parent] :=
        houtputFrame.2 parent hparentNeChild (by omega) (by omega)
      exact ⟨output, houtputRun,
        hframeRightReady.trans_fresh houtputFrame (by
          dsimp [child]
          omega),
        by omega,
        by
          dsimp [mid] at hrightCountEq houtputSizeExact ⊢
          have hleftCount :
              henselTreeInternalNodeCount start ((start + stop) / 2) = 0 := by
            rw [henselTreeInternalNodeCount, dif_neg hleft]
          simp [pushed] at hrightReadySizeEq
          omega,
        by
          have hchildEq :
              child = parent + 1 +
                henselTreeInternalNodeCount start ((start + stop) / 2) := by
            dsimp [child, mid]
            have hleftCount :
                henselTreeInternalNodeCount start ((start + stop) / 2) = 0 := by
              rw [henselTreeInternalNodeCount, dif_neg hleft]
            omega
          simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
            houtputParent, hrightReadySelected, hpushedParent,
            hafterLeftSelected, hchildEq],
        hstoredInvariant.of_algebraEq (by
          simpa [houtputParent] using hrightReadyStored),
        by
          have hrightFinal := houtputCertificate.lower_mono (by
            dsimp [child]
            omega : parent ≤ child)
          have hchildEq : child = parent + 1 +
              henselTreeInternalNodeCount start ((start + stop) / 2) := by
            dsimp [child, mid]
            have hleftCount : henselTreeInternalNodeCount start
                ((start + stop) / 2) = 0 := by
              rw [henselTreeInternalNodeCount, dif_neg hleft]
            omega
          exact .node start stop parent output output[parent] (by omega)
            (@Array.getElem?_eq_getElem _ output parent (by omega))
            (by
              simpa using (show HenselTreeNodeRawTopologyMatches start stop parent
                output[parent] from by
                  simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
                    houtputParent, hrightReadySelected, hpushedParent,
                    hafterLeftSelected, hchildEq]))
            (fun h => False.elim (hleft h))
            (fun _ => by simpa [hchildEq] using hrightFinal),
        by
          have hrightFinal := houtputSemantic.lower_mono (by
            dsimp [child]
            omega : parent ≤ child)
          have hchildEq : child = parent + 1 +
              henselTreeInternalNodeCount start ((start + stop) / 2) := by
            dsimp [child, mid]
            have hleftCount : henselTreeInternalNodeCount start
                ((start + stop) / 2) = 0 := by
              rw [henselTreeInternalNodeCount, dif_neg hleft]
            omega
          rw [henselTreeBuildTopology]
          rw [if_neg hleft, if_pos hright]
          exact .node start stop parent _ _ output output[parent] (by omega)
            (@Array.getElem?_eq_getElem _ output parent (by omega))
            (hstoredInvariant.of_algebraEq (by
              simpa [houtputParent] using hrightReadyStored))
            (hstoredOneHead.of_algebraEq (by
              simpa [houtputParent] using hrightReadyStored))
            (fun childTree hchildTree => by contradiction)
            (fun childTree hchildTree => by
              injection hchildTree with htree
              subst childTree
              simpa [hchildEq] using hrightFinal),
        houtputCanonical⟩
    · rw [dif_neg hright]
      have hparentAfterLeft : parent < afterLeft.size := by omega
      rcases henselTreeSetNodeRawIR_preserves_algebra afterLeft parent
          (fun node => { node with right := -1 }) hparentAfterLeft
          (by simp [HenselTreeNodeAlgebraEq]) with
        ⟨output, houtputRun, houtputSizeFromRun, houtputFrame,
          houtputParentBound,
          houtputAlgebra, houtputSelected⟩
      have houtputCanonical : HenselArrayCanonical output :=
        houtputFrame.canonical_of_size_eq houtputSizeFromRun
          hafterLeftCanonical houtputParentBound
          (hstoredCanonical.of_algebraEq
            (houtputAlgebra.trans hafterLeftStored))
      exact ⟨output, houtputRun, hframeAfterLeft.trans houtputFrame,
        by omega,
        by
          dsimp [mid] at *
          have hleftCount :
              henselTreeInternalNodeCount start ((start + stop) / 2) = 0 := by
            rw [henselTreeInternalNodeCount, dif_neg hleft]
          have hrightCount :
              henselTreeInternalNodeCount ((start + stop) / 2) stop = 0 := by
            rw [henselTreeInternalNodeCount, dif_neg hright]
          omega,
        by
          simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
            houtputSelected, hafterLeftSelected],
        hstoredInvariant.of_algebraEq
          (houtputAlgebra.trans hafterLeftStored),
        by
          exact .node start stop parent output output[parent] (by omega)
            (Array.getElem?_eq_getElem houtputParentBound)
            (by
              simpa using (show HenselTreeNodeRawTopologyMatches start stop parent
                output[parent] from by
                  simp [HenselTreeNodeRawTopologyMatches, hleft, hright,
                    houtputSelected, hafterLeftSelected]))
            (fun h => False.elim (hleft h))
            (fun h => False.elim (hright h)),
        by
          rw [henselTreeBuildTopology]
          rw [if_neg hleft, if_neg hright]
          exact .node start stop parent _ _ output output[parent] (by omega)
            (Array.getElem?_eq_getElem houtputParentBound)
            (hstoredInvariant.of_algebraEq
              (houtputAlgebra.trans hafterLeftStored))
            (hstoredOneHead.of_algebraEq
              (houtputAlgebra.trans hafterLeftStored))
            (fun childTree hchildTree => by contradiction)
            (fun childTree hchildTree => by contradiction),
        houtputCanonical⟩
termination_by nodes start stop parent _ _ _ => stop - start
decreasing_by all_goals omega

/-- Raw-to-safe bridge for the complete generated Hensel-tree constructor.
The root allocation and the call over the full factor interval are exactly
the generated L1 entry; no mathematical tree is supplied as an oracle. -/
theorem strictHenselTreeBuildRawIR_succeeds
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (htwo : 2 ≤ factors.size) :
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (strictHenselTreeBuildRawOps this mulProvider) factors this._p =
        .ok output ∧
      1 ≤ output.size := by
  rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
      mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
      #[default] 0 factors.size
      0 (by simpa using htwo) (by simp) (by simp) (by simp)
      HenselArrayCanonical.singleton_default with
    ⟨output, hrun, hframe, hrootBound, hsizeExact, hrootTopology,
      hrootInvariant, hrawCertificate, hsemanticCertificate, hcanonical⟩
  refine ⟨output, ?_, by simpa using hframe.1⟩
  simpa [Generated.StrictHensel.__hensel_tree_build_raw_ir, htwo] using hrun

/-- Semantic raw-to-safe entry theorem for the generated tree constructor.
In addition to executing the real root allocation and recursion, it exposes
the exact interval-product gcd certificate stored at the concrete root. -/
theorem strictHenselTreeBuildRawIR_refines_gcd
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (htwo : 2 ≤ factors.size) :
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (strictHenselTreeBuildRawOps this mulProvider) factors this._p =
        .ok output ∧
      ∃ hroot : 0 < output.size,
      HenselTreeNodeGCDInvariant this._p.toNat factors 0 factors.size
        (getElem output 0 hroot) := by
  rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
      mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
      #[default] 0 factors.size
      0 (by simpa using htwo) (by simp) (by simp) (by simp)
      HenselArrayCanonical.singleton_default with
    ⟨output, hrun, hframe, hrootBound, hsizeExact, htopology, hinvariant,
      hrawCertificate, hsemanticCertificate, hcanonical⟩
  have hroot : 0 < output.size := by simpa using hframe.1
  exact ⟨output, by
    simpa [Generated.StrictHensel.__hensel_tree_build_raw_ir, htwo] using hrun,
    hroot, hinvariant⟩

/-- The actual number of allocated nodes is exactly the node count of the
canonical source-derived finite topology.  This rules out both skipped and
extra recursive allocations. -/
theorem strictHenselTreeBuildRawIR_refines_topology_size
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (htwo : 2 ≤ factors.size) :
    let tree := henselTreeBuildTopology 0 factors.size 0
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (strictHenselTreeBuildRawOps this mulProvider) factors this._p =
        .ok output ∧
      ∃ hroot : 0 < output.size,
      output.size = tree.nodeCount ∧
      tree.rootIndex = 0 ∧
      HenselTreeNodeGCDInvariant this._p.toNat factors 0 factors.size
        (getElem output 0 hroot) := by
  dsimp only
  rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
      mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
      #[default] 0 factors.size
      0 (by simpa using htwo) (by simp) (by simp) (by simp)
      HenselArrayCanonical.singleton_default with
    ⟨output, hrun, hframe, hrootBound, hsizeExact, hrawTopology, hinvariant,
      hrawCertificate, hsemanticCertificate, hcanonical⟩
  have hcount := henselTreeInternalNodeCount_eq 0 factors.size
    (by simpa using htwo)
  have htopologyCount := henselTreeBuildTopology_nodeCount 0 factors.size 0
    (by simpa using htwo)
  have hroot : 0 < output.size := by simpa using hframe.1
  refine ⟨output, ?_, hroot, ?_, by simp, hinvariant⟩
  · simpa [Generated.StrictHensel.__hensel_tree_build_raw_ir, htwo] using hrun
  · rw [htopologyCount, hcount]
    simpa using hsizeExact

/-- With coprime root interval products, the concrete root already satisfies
the unit Bézout invariant required by the first Hensel lifting iteration. -/
theorem strictHenselTreeBuildRawIR_refines_initial_root
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (htwo : 2 ≤ factors.size)
    (hcoprime : IsCoprime
      (henselFactorRangeProduct this._p.toNat factors
        (factors.size / 2) 0)
      (henselFactorRangeProduct this._p.toNat factors factors.size
        (factors.size / 2))) :
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (strictHenselTreeBuildRawOps this mulProvider) factors this._p =
        .ok output ∧
      ∃ hroot : 0 < output.size,
      HenselTreeNodeInitialInvariant this._p.toNat factors 0 factors.size
        (getElem output 0 hroot) := by
  rcases strictHenselTreeBuildRawIR_refines_gcd this hcfg h2p hp2 mulProvider
      factors hfactors hfactorsNonempty hfactorsMonicAfterZero htwo with
    ⟨output, hrun, hroot, hinvariant⟩
  exact ⟨output, hrun, hroot,
    hinvariant.toInitial (by simpa using hcoprime)⟩

/-- The root's unit Bézout certificate follows solely from pairwise
coprimality of the concrete input factors; no per-node coprimality oracle is
required. -/
theorem strictHenselTreeBuildRawIR_refines_initial_root_of_pairwise
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (hpairwise : ∀ i j (hi : i < factors.size) (hj : j < factors.size),
      i < j → IsCoprime
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat (getElem factors i hi))
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat (getElem factors j hj)))
    (htwo : 2 ≤ factors.size) :
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (strictHenselTreeBuildRawOps this mulProvider) factors this._p =
        .ok output ∧
      ∃ hroot : 0 < output.size,
      HenselTreeNodeInitialInvariant this._p.toNat factors 0 factors.size
        (getElem output 0 hroot) := by
  apply strictHenselTreeBuildRawIR_refines_initial_root this hcfg h2p hp2
    mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero htwo
  exact henselFactorRangeProducts_isCoprime this._p.toNat factors hpairwise
    0 (factors.size / 2) factors.size (by omega) (by omega) (by simp)

/-- Root-level topology and semantic bridge for the exact generated builder.
The signed-index bound is precisely the source `h_fits_int32` requirement. -/
theorem strictHenselTreeBuildRawIR_refines_topology_root
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hfactorsMonicAfterZero : ∀ index (hindex : index < factors.size),
      0 < index → (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
        (getElem factors index hindex)).Monic)
    (hpairwise : ∀ i j (hi : i < factors.size) (hj : j < factors.size),
      i < j → IsCoprime
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat (getElem factors i hi))
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat (getElem factors j hj)))
    (htwo : 2 ≤ factors.size)
    (hfitsInt32 : henselTreeInternalNodeCount 0 factors.size < 2 ^ 31) :
    let tree := henselTreeBuildTopology 0 factors.size 0
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (strictHenselTreeBuildRawOps this mulProvider) factors this._p =
        .ok output ∧
      ∃ hroot : 0 < output.size,
      output.size = tree.nodeCount ∧
      tree.rootIndex = 0 ∧
      liftChildMatches (getElem output 0 hroot).left
          (match tree with | .node _ left _ => left) ∧
      liftChildMatches (getElem output 0 hroot).right
          (match tree with | .node _ _ right => right) ∧
      HenselExtractInvariant tree output ∧
      HenselTreeSemanticBuildCertificate this._p.toNat factors 0
        0 factors.size tree output ∧
      HenselTreeNodeInitialInvariant this._p.toNat factors 0 factors.size
        (getElem output 0 hroot) ∧
      HenselArrayCanonical output := by
  dsimp only
  rcases strictHenselTreeBuildRecursiveRawIR_succeeds this hcfg h2p hp2
      mulProvider factors hfactors hfactorsNonempty hfactorsMonicAfterZero
      #[default] 0 factors.size
      0 (by simpa using htwo) (by simp) (by simp) (by simp)
      HenselArrayCanonical.singleton_default with
    ⟨output, hrun, hframe, hrootBound, hsizeExact, hrawTopology,
      hgcdInvariant, hrawCertificate, hsemanticCertificate,
      hbuilderCanonical⟩
  have hcount := henselTreeInternalNodeCount_eq 0 factors.size
    (by simpa using htwo)
  have htopologyCount := henselTreeBuildTopology_nodeCount 0 factors.size 0
    (by simpa using htwo)
  have hchildren := hrawTopology.liftChildMatches
    (by simpa using htwo) (by simpa using hfitsInt32)
  have hextract : HenselExtractInvariant
      (henselTreeBuildTopology 0 factors.size 0) output :=
    (hrawCertificate.toExtractCertificate (by simpa using htwo) (by
      simpa using hfitsInt32)).toInvariant
  have hcoprime := henselFactorRangeProducts_isCoprime this._p.toNat factors
    hpairwise 0 (factors.size / 2) factors.size (by omega) (by omega) (by simp)
  have hroot : 0 < output.size := by simpa using hframe.1
  refine ⟨output, ?_, hroot, ?_, by simp, ?_, ?_, hextract,
    hsemanticCertificate, ?_, hbuilderCanonical⟩
  · simpa [Generated.StrictHensel.__hensel_tree_build_raw_ir, htwo] using hrun
  · rw [htopologyCount, hcount]
    simpa using hsizeExact
  · exact hchildren.1
  · exact hchildren.2
  · exact hgcdInvariant.toInitial (by simpa using hcoprime)

/-- Safety and inductive premises for the complete generated C++ Hensel
entry.  Every premise is universally quantified over the intermediate value
actually returned by the preceding raw stage; no field can prescribe a tree,
node array, extracted factor array, modulus, or final answer. -/
structure HenselLiftEntryInvariant
    (this : DenseUPolyZp)
    (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp)
    (aTarget : Int32) : Prop where
  sourceLeading : ∃ leading, f[0]? = some leading
  exponent : aTarget = 0 ∨ aTarget > 0
  targetNonnegative : ∀ target,
    HenselLiftTargetCorrect f this._p aTarget target → 0 ≤ target
  factorCount : 2 ≤ factors.size
  topologyFits : henselTreeInternalNodeCount 0 factors.size < 2 ^ 31
  adjustedCanonical : ∀ adjusted,
    HenselAdjustFirstFactorCorrect f factors this._p adjusted →
    ∀ factor ∈ adjusted.toList,
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat factor
  adjustedNonempty : ∀ adjusted,
    HenselAdjustFirstFactorCorrect f factors this._p adjusted →
    ∀ factor ∈ adjusted.toList, 0 < factor.size
  factorsMonic : ∀ index (hindex : index < factors.size),
    (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
      (getElem factors index hindex)).Monic
  adjustedPairwise : ∀ adjusted,
    HenselAdjustFirstFactorCorrect f factors this._p adjusted →
    ∀ i j (hi : i < adjusted.size) (hj : j < adjusted.size), i < j →
      IsCoprime
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
          (getElem adjusted i hi))
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat
          (getElem adjusted j hj))
  liftReady : ∀ target adjusted nodes,
    HenselLiftTargetCorrect f this._p aTarget target →
    HenselAdjustFirstFactorCorrect f factors this._p adjusted →
    Generated.StrictHensel.__hensel_tree_build_raw_ir
        (strictHenselTreeBuildRawOps this mulProvider) adjusted this._p =
          .ok nodes →
    HenselLiftLoopRefinementInvariant termination
      (henselTreeBuildTopology 0 factors.size 0) f target.toNat this._p.toNat
      nodes
  normalizeReady : ∀ target adjusted nodes liftedNodes outputM extracted,
    HenselLiftTargetCorrect f this._p aTarget target →
    HenselAdjustFirstFactorCorrect f factors this._p adjusted →
    Generated.StrictHensel.__hensel_tree_build_raw_ir
        (strictHenselTreeBuildRawOps this mulProvider) adjusted this._p =
          .ok nodes →
    HenselLiftLoopCorrect termination
      (henselTreeBuildTopology 0 factors.size 0) f target.toNat this._p.toNat
      nodes liftedNodes outputM →
    HenselExtractCorrect (henselTreeBuildTopology 0 factors.size 0)
      liftedNodes #[] extracted →
    HenselNormalizeExecutionInvariant extracted outputM

/-- Complete semantic trace of the original C++ `__hensel_lift_upoly`
entry, retaining the correctness evidence for every executed stage. -/
def HenselLiftEntryCorrect
    (termination : Generated.StrictHensel.DivmodTermination)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64)
    [Fact (Nat.Prime p.toNat)]
    (aTarget : Int32) (output : Array SparsePolyZZ × ZZ) : Prop :=
  ∃ target adjusted nodes liftedNodes outputM extracted,
    HenselLiftTargetCorrect f p aTarget target ∧
    HenselAdjustFirstFactorCorrect f factors p adjusted ∧
    HenselTreeSemanticBuildCertificate p.toNat adjusted 0
      0 factors.size (henselTreeBuildTopology 0 factors.size 0) nodes ∧
    HenselLiftLoopCorrect termination
      (henselTreeBuildTopology 0 factors.size 0) f target.toNat p.toNat nodes
        liftedNodes outputM ∧
    HenselExtractCorrect (henselTreeBuildTopology 0 factors.size 0)
      liftedNodes #[] extracted ∧
    HenselNormalizeCorrect extracted outputM output.1 ∧
    output.2 = outputM ∧
    HenselFactorArrayCanonical output.1

private theorem henselLeafOrigins_trans
    {p : Nat} {inputs : List SparsePolyZp} {initial final : List SparsePolyZZ}
    (horigins : List.Forall₂
      (fun input lifted => toPolyMod p lifted =
        CLPoly.Math.SparsePolyZp.toPoly p input) inputs initial)
    (hreduces : List.Forall₂
      (fun old new => toPolyMod p new = toPolyMod p old) initial final) :
    List.Forall₂
      (fun input lifted => toPolyMod p lifted =
        CLPoly.Math.SparsePolyZp.toPoly p input) inputs final := by
  induction horigins generalizing final with
  | nil =>
      cases hreduces
      exact .nil
  | cons horigin horigins ih =>
      cases hreduces with
      | cons hreduce hreduces =>
          exact .cons (hreduce.trans horigin) (ih hreduces)

/-- Before the source normalization block, every factor extracted by the
actual Hensel entry still reduces to its concrete adjusted finite-field input
factor, pointwise and in source order. -/
theorem HenselLiftEntryCorrect.preNormalizationOrigins
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hcount : 2 ≤ factors.size)
    (hcorrect : HenselLiftEntryCorrect termination f factors p aTarget output) :
    ∃ adjusted extracted outputM,
      HenselAdjustFirstFactorCorrect f factors p adjusted ∧
      HenselNormalizeCorrect extracted outputM output.1 ∧
      output.2 = outputM ∧
      List.Forall₂
        (fun input lifted => toPolyMod p.toNat lifted =
          CLPoly.Math.SparsePolyZp.toPoly p.toNat input)
        (henselFactorRangeList adjusted factors.size 0) extracted.toList ∧
      HenselNormalizeUnitRel p.toNat extracted output.1 := by
  rcases hcorrect with
    ⟨target, adjusted, nodes, liftedNodes, outputM, extracted,
      htarget, hadjust, hsemantic, hlift, hextract, hnormalize, houtputM,
      _houtputCanonical⟩
  have hinitial := hsemantic.extractedFactors_forall₂ hcount
  have harray := hlift.arrayReduces_of_dvd (dvd_refl p.toNat)
  have hlifted := hsemantic.extractedFactors_forall₂_of_reduces harray
  have horigins := henselLeafOrigins_trans hinitial hlifted
  have hnormalizeRel := hnormalize.unitRel hlift.initialM_dvd_outputM
  have hextracted : extracted.toList =
      henselExtractedFactors
        (henselTreeBuildTopology 0 factors.size 0) liftedNodes := by
    simpa using hextract.toList_eq
  rw [← hextracted] at horigins
  exact ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
    horigins, hnormalizeRel⟩

/-- Before the source normalization block, the factors extracted by the
actual generated Hensel entry multiply to the input polynomial at the exact
modulus returned by the well-founded lift loop.  The finite-field source
equation remains an explicit premise and is supplied by the preceding Zp
factorization refinement; Hensel lifting does not assume a semantic oracle. -/
theorem HenselLiftEntryCorrect.preNormalizationProduct
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hcount : 2 ≤ factors.size)
    (hsource : ∀ adjusted,
      HenselAdjustFirstFactorCorrect f factors p adjusted →
      ((henselFactorRangeList adjusted factors.size 0).map
        (CLPoly.Math.SparsePolyZp.toPoly p.toNat)).prod =
          toPolyMod p.toNat f)
    (hcorrect : HenselLiftEntryCorrect termination f factors p aTarget output) :
    ∃ adjusted, ∃ extracted, ∃ outputM : Nat,
      HenselAdjustFirstFactorCorrect f factors p adjusted ∧
      HenselNormalizeCorrect extracted (outputM : Int) output.1 ∧
      output.2 = (outputM : Int) ∧
      (extracted.toList.map (toPolyMod outputM)).prod =
        toPolyMod outputM f := by
  rcases hcorrect with
    ⟨target, adjusted, nodes, liftedNodes, outputM, extracted,
      _htarget, hadjust, hsemantic, hlift, hextract, hnormalize, houtputM,
      _houtputCanonical⟩
  have hproduct := hsemantic.liftLoop_extractedFactors_product hlift hcount
    (hsource adjusted hadjust)
  have hextracted : extracted.toList =
      henselExtractedFactors
        (henselTreeBuildTopology 0 factors.size 0) liftedNodes := by
    simpa using hextract.toList_eq
  have hproduct' :
      (extracted.toList.map (toPolyMod outputM)).prod =
        toPolyMod outputM f := by
    rw [hextracted]
    exact hproduct
  exact ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
    hproduct'⟩

/-- A successful full generated Hensel entry exposes the precision reached by
its own well-founded loop.  In particular, the returned modulus is not an
oracle-provided value: it strictly exceeds the concrete target computed by
`__hensel_lift_target_raw_ir`. -/
theorem HenselLiftEntryCorrect.outputModulus_gt_target
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (htargetNonnegative : ∀ target,
      HenselLiftTargetCorrect f p aTarget target → 0 ≤ target)
    (hcorrect : HenselLiftEntryCorrect termination f factors p aTarget output) :
    ∃ target : ZZ,
      HenselLiftTargetCorrect f p aTarget target ∧ target < output.2 := by
  rcases hcorrect with
    ⟨target, adjusted, nodes, liftedNodes, outputM, extracted,
      htarget, hadjust, _hsemantic, hlift, hextract, hnormalize, houtputM,
      _houtputCanonical⟩
  refine ⟨target, htarget, ?_⟩
  rw [houtputM]
  have htargetNonnegative' := htargetNonnegative target htarget
  have hprecision := hlift.outputM_gt_target
  calc
    target = (target.toNat : Int) := (Int.toNat_of_nonneg htargetNonnegative').symm
    _ < (outputM : Int) := by exact_mod_cast hprecision

private theorem henselLiftUpolyRawIR_run_of_stages
    (stepOps : Generated.StrictHensel.HenselStepRawOps)
    (treeOps : Generated.StrictHensel.HenselTreeBuildRawOps)
    (f : SparsePolyZZ) (factors adjusted : Array SparsePolyZp)
    (p : UInt64) (aTarget : Int32) (hp : 2 ≤ p.toNat)
    (target : ZZ) (nodes liftedNodes : Array HenselNode)
    (outputM : Nat) (extracted normalized : Array SparsePolyZZ)
    (hcount : 2 ≤ factors.size)
    (htarget : Generated.StrictHensel.__hensel_lift_target_raw_ir
      f p aTarget = .ok target)
    (htargetNonnegative : 0 ≤ target)
    (hadjust : Generated.StrictHensel.__hensel_adjust_first_factor_raw_ir
      f factors p = .ok adjusted)
    (hbuild : Generated.StrictHensel.__hensel_tree_build_raw_ir
      treeOps adjusted p = .ok nodes)
    (hlift : Generated.StrictHensel.__hensel_lift_loop_raw_ir stepOps
      (Generated.StrictHensel.henselTreeBuildTopologyRawIR 0 factors.size 0)
      f target.toNat p.toNat hp nodes = .ok (liftedNodes, outputM))
    (hextract : Generated.StrictHensel.__hensel_extract_factors_raw_ir
      (Generated.StrictHensel.henselTreeBuildTopologyRawIR 0 factors.size 0)
      liftedNodes #[] = .ok extracted)
    (hnormalize : Generated.StrictHensel.__hensel_normalize_result_raw_ir
      extracted outputM = .ok normalized) :
    Generated.StrictHensel.__hensel_lift_upoly_raw_ir stepOps treeOps
      f factors p aTarget hp = .ok (normalized, outputM) := by
  rw [Generated.StrictHensel.__hensel_lift_upoly_raw_ir]
  simp only [hcount, ↓reduceIte, pure, Except.pure, htarget, bind,
    Except.bind]
  split
  · simp only [hadjust, Except.bind]
    rw [hbuild]
    simp only [Except.bind]
    rw [hlift]
    simp only [Except.bind]
    rw [hextract]
    simp only [Except.bind]
    rw [hnormalize]
  · contradiction

/-- Genuine raw-to-safe and L1-to-L2 composition theorem for the full C++
Hensel entry.  Every intermediate is obtained from the strict generated raw
program before its invariant is instantiated. -/
theorem __hensel_lift_upoly_raw_ir_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (aTarget : Int32)
    (hinvariant : HenselLiftEntryInvariant this termination mulProvider f
      factors aTarget) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
          (strictHenselRawOps termination)
          (strictHenselTreeBuildRawOps this mulProvider)
          f factors this._p aTarget
            (Nat.Prime.two_le (Fact.out : Nat.Prime this._p.toNat)) =
              .ok output ∧
      HenselLiftEntryCorrect termination f factors this._p aTarget output := by
  rcases __hensel_lift_target_raw_ir_refines f this._p aTarget
      hinvariant.sourceLeading hinvariant.exponent with
    ⟨target, htargetRun, htargetCorrect⟩
  have htargetNonnegative := hinvariant.targetNonnegative target htargetCorrect
  rcases __hensel_adjust_first_factor_raw_ir_refines f factors this._p
      ⟨hinvariant.sourceLeading, by
        have hzero : 0 < factors.size := lt_of_lt_of_le (by omega)
          hinvariant.factorCount
        exact ⟨factors[0], Array.getElem?_eq_getElem hzero⟩⟩ with
    ⟨adjusted, hadjustRun, hadjustCorrect⟩
  have hadjustSize : adjusted.size = factors.size := by
    cases hadjustCorrect with
    | adjusted leading first value hsource hfirst hvalue =>
        have hindex : 0 < factors.size := by
          by_contra hnot
          rw [Array.getElem?_eq_none (by omega)] at hfirst
          contradiction
        simp [Array.set!, hindex]
  let tree := henselTreeBuildTopology 0 factors.size 0
  rcases strictHenselTreeBuildRawIR_refines_topology_root this hcfg h2p hp2
      mulProvider adjusted (hinvariant.adjustedCanonical adjusted hadjustCorrect)
      (hinvariant.adjustedNonempty adjusted hadjustCorrect)
      (hadjustCorrect.monic_of_pos hinvariant.factorsMonic)
      (hinvariant.adjustedPairwise adjusted hadjustCorrect)
      (by rw [hadjustSize]; exact hinvariant.factorCount)
      (by rw [hadjustSize]; exact hinvariant.topologyFits) with
    ⟨nodes, hnodesRun, hroot, hsize, hrootIndex, hleft, hright,
      hextractInvariant, hsemanticInvariant, hrootInvariant,
      hnodesCanonical⟩
  have htreeEq : Generated.StrictHensel.henselTreeBuildTopologyRawIR
      0 factors.size 0 = tree := henselTreeBuildTopologyRawIR_eq 0 factors.size 0
  have hsemanticInvariant' :
      HenselTreeSemanticBuildCertificate this._p.toNat adjusted 0
        0 factors.size tree nodes := by
    simpa [tree, hadjustSize] using hsemanticInvariant
  rcases __hensel_lift_loop_raw_ir_refines termination tree f target.toNat
      this._p.toNat nodes
      (hinvariant.liftReady target adjusted nodes htargetCorrect hadjustCorrect
        hnodesRun) with
    ⟨liftedNodes, outputM, hliftRun, hliftCorrect⟩
  have hextractInvariant' : HenselExtractInvariant tree nodes := by
    simpa [tree, hadjustSize] using hextractInvariant
  have hliftedExtract :=
    hextractInvariant'.of_topologyEq hliftCorrect.topologyEq
  rcases __hensel_extract_factors_raw_ir_refines tree liftedNodes #[]
      hliftedExtract with ⟨extracted, hextractRun, hextractCorrect⟩
  have hliftedCanonical : HenselArrayCanonical liftedNodes :=
    hliftCorrect.arrayCanonical hnodesCanonical
  have hemptyCanonical : HenselFactorArrayCanonical (#[] : Array SparsePolyZZ) := by
    intro factor hfactor
    simp at hfactor
  have hextractedCanonical : HenselFactorArrayCanonical extracted :=
    hextractCorrect.outputCanonical hliftedCanonical hemptyCanonical
  rcases __hensel_normalize_result_raw_ir_refines extracted outputM
      (hinvariant.normalizeReady target adjusted nodes liftedNodes outputM
        extracted htargetCorrect hadjustCorrect
        hnodesRun hliftCorrect hextractCorrect) with
    ⟨normalized, hnormalizeRun, hnormalizeCorrect⟩
  have hnormalizedCanonical : HenselFactorArrayCanonical normalized :=
    hnormalizeCorrect.outputCanonical hextractedCanonical
  refine ⟨(normalized, outputM), ?_,
    ⟨target, adjusted, nodes, liftedNodes, outputM, extracted,
      htargetCorrect, hadjustCorrect, hsemanticInvariant', hliftCorrect,
      hextractCorrect,
      hnormalizeCorrect, rfl, hnormalizedCanonical⟩⟩
  apply henselLiftUpolyRawIR_run_of_stages
    (hcount := hinvariant.factorCount) (htarget := htargetRun)
    (htargetNonnegative := htargetNonnegative) (hadjust := hadjustRun)
    (hbuild := hnodesRun)
  · simpa [htreeEq] using hliftRun
  · simpa [htreeEq] using hextractRun
  · exact hnormalizeRun

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
