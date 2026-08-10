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
      toPolyMod m factorNode.h = toPolyMod m node.h := by
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
  · simpa [factorNode] using hsemantic

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

end StrictHensel

-- No Hensel L1→L2 theorem or legacy candidate is exported until a strict
-- cpp2lean-generated entry and its direct execution proof are available.

end Refinement
