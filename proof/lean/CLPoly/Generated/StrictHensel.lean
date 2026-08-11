-- Auto-generated strict raw C++ control flow for `__hensel_step`.
import CLPoly.Model

set_option autoImplicit false
set_option maxHeartbeats 1600000

namespace Generated.StrictHensel

def certifyRawExec {α : Type} (run : RawExec α) :
    RawExec { output : α // run = .ok output } :=
  match hrun : run with
  | .error fault => .error fault
  | .ok output => .ok ⟨output, by simpa using hrun⟩

/-- Exact lowering of the source coefficient loop implementing
`coefficient := fdiv_q(coefficient,m); coefficient := fdiv_r(coefficient,m)`
followed by the in-place zero compaction. -/
def divideCoeffs (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.map fun term => (term.fst, ZZ.fdiv_q term.snd term.snd m)

def divideThenReduceCoeffs (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.filterMap fun term =>
    let quotient := ZZ.fdiv_q term.snd term.snd m
    let coefficient := ZZ.fdiv_r quotient quotient m
    if coefficient != 0 then some (term.fst, coefficient) else none

/-- Exact total lowering shared by the four source range-for loops which
multiply every coefficient by `m`. -/
def scaleCoeffs (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.map fun term => (term.fst, term.snd * m)

/-- Exact total lowering of C++ `__upoly_mod_coeff`: floor remainder on every
coefficient followed by removal of zero terms. -/
def modCoeffOutput (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.filterMap fun term =>
    let coefficient := ZZ.fdiv_r term.snd term.snd m
    if coefficient != 0 then some (term.fst, coefficient) else none

def __upoly_mod_coeff_raw_ir (f : SparsePolyZZ) (m : ZZ) :
    RawExec SparsePolyZZ :=
  .ok (modCoeffOutput f m)

def pushNonzero (result : SparsePolyZZ) (degree : Nat)
    (coefficient : ZZ) : SparsePolyZZ :=
  if coefficient != 0 then result.push (UMonomial.mk degree, coefficient)
  else result

/-- Exact lowering of `pair_vec_add` for the canonical univariate vectors
used by `basic_polynomial::operator+`.  The two source iterator positions are
explicit and every recursive branch advances at least one iterator. -/
def pairVecAddLoop (a b : SparsePolyZZ) (aIndex bIndex : Nat)
    (result : SparsePolyZZ) : SparsePolyZZ :=
  if hmore : aIndex < a.size ∨ bIndex < b.size then
    if haDone : aIndex ≥ a.size then
      pairVecAddLoop a b aIndex (bIndex + 1) (result.push b[bIndex]!)
    else if hbDone : bIndex ≥ b.size then
      pairVecAddLoop a b (aIndex + 1) bIndex (result.push a[aIndex]!)
    else if hdegree : b[bIndex]!.1.deg > a[aIndex]!.1.deg then
      pairVecAddLoop a b aIndex (bIndex + 1) (result.push b[bIndex]!)
    else if hequal : b[bIndex]!.1.deg = a[aIndex]!.1.deg then
      pairVecAddLoop a b (aIndex + 1) (bIndex + 1)
        (pushNonzero result a[aIndex]!.1.deg
          (b[bIndex]!.2 + a[aIndex]!.2))
    else
      pairVecAddLoop a b (aIndex + 1) bIndex (result.push a[aIndex]!)
  else result
termination_by (a.size - aIndex) + (b.size - bIndex)
decreasing_by all_goals simp_wf; omega

def __upoly_add_raw_ir (a b : SparsePolyZZ) : RawExec SparsePolyZZ :=
  .ok (pairVecAddLoop a b 0 0 #[])

/-- Exact lowering of `pair_vec_sub`.  Terms consumed only from the right
source are negated, matching the C++ implementation. -/
def pairVecSubLoop (a b : SparsePolyZZ) (aIndex bIndex : Nat)
    (result : SparsePolyZZ) : SparsePolyZZ :=
  if hmore : aIndex < a.size ∨ bIndex < b.size then
    if haDone : aIndex ≥ a.size then
      pairVecSubLoop a b aIndex (bIndex + 1)
        (result.push (b[bIndex]!.1, -b[bIndex]!.2))
    else if hbDone : bIndex ≥ b.size then
      pairVecSubLoop a b (aIndex + 1) bIndex (result.push a[aIndex]!)
    else if hdegree : b[bIndex]!.1.deg > a[aIndex]!.1.deg then
      pairVecSubLoop a b aIndex (bIndex + 1)
        (result.push (b[bIndex]!.1, -b[bIndex]!.2))
    else if hequal : b[bIndex]!.1.deg = a[aIndex]!.1.deg then
      pairVecSubLoop a b (aIndex + 1) (bIndex + 1)
        (pushNonzero result a[aIndex]!.1.deg
          (a[aIndex]!.2 - b[bIndex]!.2))
    else
      pairVecSubLoop a b (aIndex + 1) bIndex (result.push a[aIndex]!)
  else result
termination_by (a.size - aIndex) + (b.size - bIndex)
decreasing_by all_goals simp_wf; omega

def __upoly_sub_raw_ir (a b : SparsePolyZZ) : RawExec SparsePolyZZ :=
  .ok (pairVecSubLoop a b 0 0 #[])

/-- One pending `addmul` contribution in the C++ multiplication heap. -/
structure MulProduct where
  degree : Nat
  coefficient : ZZ
deriving DecidableEq

def pairVecMulProducts (a b : SparsePolyZZ) : List MulProduct :=
  a.toList.flatMap fun ta =>
    b.toList.map fun tb =>
      ⟨ta.1.deg + tb.1.deg, ta.2 * tb.2⟩

def mulMaxDegree (products : List MulProduct) : Nat :=
  (products.map (fun product => product.degree)).max?.getD 0

def mulDegreeCoefficient (degree : Nat) : List MulProduct → ZZ
  | [] => 0
  | product :: products =>
      if product.degree = degree then
        product.coefficient + mulDegreeCoefficient degree products
      else mulDegreeCoefficient degree products

/-- Semantic lowering of the observable C++ multiplication heap loop.  A
step extracts the maximum monomial, performs every `addmul` attached to that
heap key, conditionally appends the coefficient, and leaves the other heap
frontier entries for the next extraction. -/
def pairVecMulHeapLoop (products : List MulProduct)
    (result : SparsePolyZZ) : SparsePolyZZ :=
  if hempty : products = [] then result
  else
    let degree := mulMaxDegree products
    let coefficient := mulDegreeCoefficient degree products
    let remaining := products.filter fun product => product.degree != degree
    pairVecMulHeapLoop remaining (pushNonzero result degree coefficient)
termination_by products.length
decreasing_by
  simp_wf
  have hsome :
      (products.map (fun product => product.degree)).max? =
        some (mulMaxDegree products) := by
    unfold mulMaxDegree
    cases hmax : (products.map (fun product => product.degree)).max? with
    | none =>
        have := List.max?_eq_none_iff.mp hmax
        simp at this
        contradiction
    | some maximum => simp [hmax]
  have hdegreeMem : mulMaxDegree products ∈
      products.map (fun product => product.degree) :=
    List.max?_mem hsome
  have hexists : ∃ product ∈ products,
      product.degree = mulMaxDegree products := by
    exact List.mem_map.mp hdegreeMem
  exact hexists

def __upoly_mul_raw_ir (a b : SparsePolyZZ) : RawExec SparsePolyZZ :=
  .ok (pairVecMulHeapLoop (pairVecMulProducts a b) #[])

/-- The inner ordered merge implementing
`r -= coefficient * x^degreeShift * g`.  Each branch advances at least one
source iterator, so termination is structural on the two remaining suffixes. -/
def divmodMergeLoop (r g : SparsePolyZZ) (coefficient m : ZZ)
    (degreeShift rIndex gIndex : Nat) (result : SparsePolyZZ) :
    SparsePolyZZ :=
  if hmore : rIndex < r.size ∨ gIndex < g.size then
    if hgDone : gIndex ≥ g.size then
      let c := ZZ.fdiv_r r[rIndex]!.2 r[rIndex]!.2 m
      divmodMergeLoop r g coefficient m degreeShift (rIndex + 1) gIndex
        (pushNonzero result r[rIndex]!.1.deg c)
    else if hrDone : rIndex ≥ r.size then
      let degree := g[gIndex]!.1.deg + degreeShift
      let cRaw := m - ((coefficient * g[gIndex]!.2 % m) % m)
      let c := ZZ.fdiv_r cRaw cRaw m
      divmodMergeLoop r g coefficient m degreeShift rIndex (gIndex + 1)
        (pushNonzero result degree c)
    else
      let rDegree := r[rIndex]!.1.deg
      let gDegree := g[gIndex]!.1.deg + degreeShift
      if rDegree > gDegree then
        let c := ZZ.fdiv_r r[rIndex]!.2 r[rIndex]!.2 m
        divmodMergeLoop r g coefficient m degreeShift (rIndex + 1) gIndex
          (pushNonzero result rDegree c)
      else if rDegree < gDegree then
        let cRaw := m - ((coefficient * g[gIndex]!.2 % m) % m)
        let c := ZZ.fdiv_r cRaw cRaw m
        divmodMergeLoop r g coefficient m degreeShift rIndex (gIndex + 1)
          (pushNonzero result gDegree c)
      else
        let cRaw := r[rIndex]!.2 - coefficient * g[gIndex]!.2
        let c := ZZ.fdiv_r cRaw cRaw m
        divmodMergeLoop r g coefficient m degreeShift (rIndex + 1)
          (gIndex + 1) (pushNonzero result rDegree c)
  else result
termination_by (r.size - rIndex) + (g.size - gIndex)
decreasing_by all_goals simp_wf; omega

def divmodCoefficient (r : SparsePolyZZ) (inverse m : ZZ) : ZZ :=
  ZZ.fdiv_r (r[0]!.2 * inverse) (r[0]!.2 * inverse) m

def divmodRemainder (r g : SparsePolyZZ) (coefficient m : ZZ)
    (degreeShift : Nat) : SparsePolyZZ :=
  divmodMergeLoop r g coefficient m degreeShift 1 1 #[]

/-- Exact finite source trace for the deterministic modular long-division
while-loop.  Each constructor records one concrete source branch. -/
inductive DivmodTrace (g : SparsePolyZZ) (inverse m : ZZ) :
    SparsePolyZZ → SparsePolyZZ → Type
  | done (r q) (inactive : ¬(!r.isEmpty && get_deg r ≥ get_deg g)) :
      DivmodTrace g inverse m r q
  | vanished (r q) (active : !r.isEmpty && get_deg r ≥ get_deg g)
      (zero : divmodCoefficient r inverse m = 0)
      (next : DivmodTrace g inverse m (r.eraseIdxIfInBounds 0) q) :
      DivmodTrace g inverse m r q
  | subtract (r q) (active : !r.isEmpty && get_deg r ≥ get_deg g)
      (nonzero : divmodCoefficient r inverse m ≠ 0)
      (next : DivmodTrace g inverse m
        (divmodRemainder r g (divmodCoefficient r inverse m) m
          (get_deg r - get_deg g).toNatClampNeg)
        (q.push (UMonomial.mk (get_deg r - get_deg g).toNatClampNeg,
          divmodCoefficient r inverse m))) :
      DivmodTrace g inverse m r q

def divmodLoop (g : SparsePolyZZ) (inverse m : ZZ) {r q : SparsePolyZZ}
    (trace : DivmodTrace g inverse m r q) : SparsePolyZZ × SparsePolyZZ :=
  match trace with
  | .done r q _ => (q, r)
  | .vanished _ _ _ _ next => divmodLoop g inverse m next
  | .subtract _ _ _ _ next => divmodLoop g inverse m next

structure DivmodTermination where
  trace : ∀ (f g reduced : SparsePolyZZ) (m : ZZ),
    __upoly_mod_coeff_raw_ir f m = .ok reduced →
    let inverse := (ZZ.invert 0 g[0]!.2 m).2
    DivmodTrace g inverse m reduced #[]

/-- Strict generated C++ `__upoly_divmod_mod` entry.  Failed source asserts
remain explicit; a successful path consumes only its exact finite trace. -/
def __upoly_divmod_mod_raw_ir (termination : DivmodTermination)
    (f g : SparsePolyZZ) (m : ZZ) : RawExec (SparsePolyZZ × SparsePolyZZ) :=
  if g.isEmpty then .error .assertionFailure
  else
    match hreduce : __upoly_mod_coeff_raw_ir f m with
    | .error fault => .error fault
    | .ok reduced =>
      let inverseRun := ZZ.invert 0 g[0]!.2 m
      if inverseRun.1 then
        .ok (divmodLoop g inverseRun.2 m
          (termination.trace f g reduced m hreduce))
      else .error .assertionFailure

structure HenselStepRawOps where
  mul : SparsePolyZZ → SparsePolyZZ → RawExec SparsePolyZZ
  add : SparsePolyZZ → SparsePolyZZ → RawExec SparsePolyZZ
  sub : SparsePolyZZ → SparsePolyZZ → RawExec SparsePolyZZ
  divmodTermination : DivmodTermination

/-- First contiguous source phase of `__hensel_step`: correct `g` and `h`
from the factorization error, while leaving the Bezout coefficients intact. -/
def __hensel_step_factor_phase_raw_ir (ops : HenselStepRawOps)
    (node : HenselNode)
    (f : SparsePolyZZ) (m : ZZ) : RawExec HenselNode := do
  let m2 := m * m
  let gh ← ops.mul node.g node.h
  let difference ← ops.sub f gh
  let e := divideThenReduceCoeffs difference m
  let se ← ops.mul node.s e
  let qr ← __upoly_divmod_mod_raw_ir ops.divmodTermination se node.h m
  let te ← ops.mul node.t e
  let qg ← ops.mul qr.1 node.g
  let tauRaw ← ops.add te qg
  let tau ← __upoly_mod_coeff_raw_ir tauRaw m
  let gRaw ← ops.add node.g (scaleCoeffs tau m)
  let gNew ← __upoly_mod_coeff_raw_ir gRaw m2
  let hRaw ← ops.add node.h (scaleCoeffs qr.2 m)
  let hNew ← __upoly_mod_coeff_raw_ir hRaw m2
  return { node with g := gNew, h := hNew }

/-- Second contiguous source phase of `__hensel_step`: correct `s` and `t`
for the already-updated factor node. -/
def __hensel_step_bezout_phase_raw_ir (ops : HenselStepRawOps)
    (factorNode : HenselNode) (m : ZZ) : RawExec HenselNode := do
  let m2 := m * m
  let sg ← ops.mul factorNode.s factorNode.g
  let th ← ops.mul factorNode.t factorNode.h
  let one : SparsePolyZZ := #[(UMonomial.mk 0, 1)]
  let oneMinusSg ← ops.sub one sg
  let bezoutDifference ← ops.sub oneMinusSg th
  let ep := divideThenReduceCoeffs bezoutDifference m
  let sep ← ops.mul factorNode.s ep
  let qrBezout ← __upoly_divmod_mod_raw_ir ops.divmodTermination
    sep factorNode.h m
  let sRaw ← ops.add factorNode.s (scaleCoeffs qrBezout.2 m)
  let sNew ← __upoly_mod_coeff_raw_ir sRaw m2
  let tep ← ops.mul factorNode.t ep
  let qpg ← ops.mul qrBezout.1 factorNode.g
  let tau2Raw ← ops.add tep qpg
  let tau2 ← __upoly_mod_coeff_raw_ir tau2Raw m
  let tRaw ← ops.add factorNode.t (scaleCoeffs tau2 m)
  let tNew ← __upoly_mod_coeff_raw_ir tRaw m2
  return { factorNode with s := sNew, t := tNew }

/-- Strict sequential translation of C++ `__hensel_step`.  It neither calls
the L2 Hensel model nor manufactures an output when a raw operation fails.
The two binds above are only proof-visible boundaries between contiguous
source statement ranges; they preserve the exact source operation order. -/
def __hensel_step_raw_ir (ops : HenselStepRawOps) (node : HenselNode)
    (f : SparsePolyZZ) (m : ZZ) : RawExec HenselNode := do
  let factorNode ← __hensel_step_factor_phase_raw_ir ops node f m
  __hensel_step_bezout_phase_raw_ir ops factorNode m

/-- Finite topology certificate for the source Hensel tree traversal.  It
contains only concrete array indices and branching shape: in particular it
contains neither polynomial values nor expected outputs.  Structural
recursion on this certificate replaces the partial recursion emitted by the
legacy translator without introducing a fuel counter. -/
inductive HenselLiftTree where
  | node (index : Nat) (left right : Option HenselLiftTree)
deriving Inhabited

def HenselLiftTree.rootIndex : HenselLiftTree → Nat
  | .node index _ _ => index

noncomputable def HenselLiftTree.nodeCount : HenselLiftTree → Nat :=
  HenselLiftTree.rec (fun _ _ _ leftCount rightCount =>
    1 + leftCount + rightCount) 0 (fun _ childCount => childCount)

theorem HenselLiftTree.left_nodeCount_lt (index : Nat)
    (child : HenselLiftTree) (right : Option HenselLiftTree) :
    child.nodeCount <
      (HenselLiftTree.node index (some child) right).nodeCount := by
  simp [HenselLiftTree.nodeCount]
  omega

theorem HenselLiftTree.right_nodeCount_lt (index : Nat)
    (left : Option HenselLiftTree) (child : HenselLiftTree) :
    child.nodeCount <
      (HenselLiftTree.node index left (some child)).nodeCount := by
  simp [HenselLiftTree.nodeCount]
  omega

/-- Exact strict lowering of C++ `__hensel_lift_recursive` on a certified
finite tree.  The program executes and stores the concrete Hensel step before
visiting the left child.  After the left traversal it re-reads the parent
node before selecting the right target, matching the two source `nodes[idx]`
reads rather than retaining a proof-side copy.

Out-of-bounds accesses remain explicit raw faults.  Agreement between the
certificate and the source `left/right == -1` branches is a refinement
precondition proved outside this generated program. -/
def __hensel_lift_recursive_raw_ir (ops : HenselStepRawOps) :
    HenselLiftTree → Array HenselNode → SparsePolyZZ → ZZ →
      RawExec (Array HenselNode)
  | .node index left right, nodes, target, m => do
      match nodes[index]? with
      | none => .error (.outOfBounds 0 index)
      | some node =>
          let lifted ← __hensel_step_raw_ir ops node target m
          let nodes := nodes.set! index lifted
          let nodes ← match left with
            | none => .ok nodes
            | some child =>
                __hensel_lift_recursive_raw_ir ops child nodes lifted.g m
          match nodes[index]? with
          | none => .error (.outOfBounds 0 index)
          | some parent =>
              match right with
              | none => .ok nodes
              | some child =>
                  __hensel_lift_recursive_raw_ir ops child nodes parent.h m

/-- Exact lowering of the quadratic-precision `while (m <= target)` loop in
C++ `__hensel_lift`.  The natural-number parameters are the nonnegative view
of the source `ZZ` values; `hm` is erased termination evidence, not fuel.
Each successful iteration executes the strict tree traversal before replacing
`m` by `m * m`, exactly as in the source. -/
def __hensel_lift_loop_raw_ir (ops : HenselStepRawOps)
    (tree : HenselLiftTree) (f : SparsePolyZZ) (target : Nat) :
    (m : Nat) → 2 ≤ m → Array HenselNode →
      RawExec (Array HenselNode × Nat)
  | m, hm, nodes =>
      if hcontinue : m ≤ target then do
        let nodes ← __hensel_lift_recursive_raw_ir ops tree nodes f (m : Int)
        __hensel_lift_loop_raw_ir ops tree f target (m * m) (by
          have hmul := Nat.mul_le_mul_left m hm
          omega) nodes
      else
        .ok (nodes, m)
termination_by m hm nodes => target + 1 - m
decreasing_by
  have hmul := Nat.mul_le_mul_left m hm
  have hgrow : m < m * m := by omega
  omega

/-- Exact strict lowering of C++ `__hensel_extract_factors`.  A missing child
pushes the corresponding concrete `g`/`h`; a present child is traversed before
the right branch, preserving source order. -/
def __hensel_extract_factors_raw_ir :
    HenselLiftTree → Array HenselNode → Array SparsePolyZZ →
      RawExec (Array SparsePolyZZ)
  | .node index left right, nodes, factors => do
      match nodes[index]? with
      | none => .error (.outOfBounds 0 index)
      | some node =>
          let factors ← match left with
            | none => .ok (factors.push node.g)
            | some child =>
                __hensel_extract_factors_raw_ir child nodes factors
          match right with
          | none => .ok (factors.push node.h)
          | some child =>
              __hensel_extract_factors_raw_ir child nodes factors

/-- Exact lowering of the final leading-coefficient normalization block in
C++ `__hensel_lift`.  Empty output takes the short-circuit branch.  A
nonempty output with an empty first polynomial exposes the source `front()`
precondition as a raw bounds fault, and a failed modular inverse exposes the
source assertion. -/
def __hensel_normalize_result_raw_ir
    (result : Array SparsePolyZZ) (m : ZZ) :
    RawExec (Array SparsePolyZZ) :=
  match result[0]? with
  | none => .ok result
  | some first =>
      match first[0]? with
      | none => .error (.outOfBounds 0 0)
      | some leading =>
          if leading.2 != 1 then
            let inverseRun := ZZ.invert 0 leading.2 m
            if inverseRun.1 then do
              let scaled := scaleCoeffs first inverseRun.2
              let normalized ← __upoly_mod_coeff_raw_ir scaled m
              .ok (result.set! 0 normalized)
            else
              .error .assertionFailure
          else
            .ok result

/-- Exact lowering of the leading-coefficient baking block at the start of
C++ `__hensel_lift`: read `lc(f)`, multiply every coefficient of factor zero
by its residue modulo `p`, normalize that factor, and write it back. -/
def scaleZpCoeffs (f : SparsePolyZp) (coefficient : Zp) : SparsePolyZp :=
  f.map fun term => (term.1, term.2 * coefficient)

def __hensel_adjust_first_factor_raw_ir
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64) :
    RawExec (Array SparsePolyZp) :=
  match f[0]? with
  | none => .error (.outOfBounds 0 0)
  | some leading =>
      match factors[0]? with
      | none => .error (.outOfBounds 1 0)
      | some first =>
          let lcModP := Zp.ofInt leading.2 p
          let adjusted := SparsePolyZp.normalization (scaleZpCoeffs first lcModP)
          .ok (factors.set! 0 adjusted)

/-- Exact finite lowering of the positive-`a_target` source loop.  `remaining`
is the number of source iterations still required, not an execution fuel:
each constructor corresponds to one `target *= p; ++i` iteration. -/
def henselExplicitTargetLoop (p : UInt64) : Nat → ZZ → ZZ
  | 0, target => target
  | remaining + 1, target =>
      henselExplicitTargetLoop p remaining (target * (p.toNat : ZZ))

/-- The `a_target > 0` target branch of C++ `__hensel_lift`, including the
final subtraction which turns `m < p^a` into `m ≤ target`. -/
def __hensel_explicit_target_raw_ir (p : UInt64) (aTarget : Int32) :
    RawExec ZZ :=
  if aTarget > 0 then
    .ok (henselExplicitTargetLoop p aTarget.toNatClampNeg 1 - 1)
  else
    .error .assertionFailure

/-- Concrete mutable state of the source extended-Euclidean loop used by
`__hensel_tree_build_recursive`. -/
structure HenselEEAState where
  r0 : SparsePolyZp
  r1 : SparsePolyZp
  s0 : SparsePolyZp
  s1 : SparsePolyZp
  t0 : SparsePolyZp
  t1 : SparsePolyZp

/-- The one still-unexpanded C++ operation at this boundary is `pair_vec_div`
with both quotient and remainder outputs.  This interface carries an
executable raw call only; it contains no expected EEA result or L2 fact. -/
structure HenselEEARawOps where
  divmod : SparsePolyZp → SparsePolyZp →
    RawExec (SparsePolyZp × SparsePolyZp)
  inverse : Zp → RawExec Zp

def henselEEAScaleNormalize (coefficient : Zp)
    (f : SparsePolyZp) : SparsePolyZp :=
  SparsePolyZp.normalization (scaleZpCoeffs f coefficient)

def henselEEANextState (state : HenselEEAState)
    (quotient remainder : SparsePolyZp) : HenselEEAState :=
  let s2 := SparsePolyZp.normalization
    (SparsePolyZp.subImpl state.s0
      (SparsePolyZp.mulImpl quotient state.s1))
  let t2 := SparsePolyZp.normalization
    (SparsePolyZp.subImpl state.t0
      (SparsePolyZp.mulImpl quotient state.t1))
  { r0 := state.r1, r1 := remainder
    s0 := state.s1, s1 := s2
    t0 := state.t1, t1 := t2 }

/-- Well-founded evidence for the actual EEA state transition.  The decrease
must hold for every successful result of the concrete raw division call, so
it cannot select a convenient quotient or remainder. -/
structure HenselEEATermination (ops : HenselEEARawOps) where
  measure : HenselEEAState → Nat
  decreases : ∀ state quotient remainder,
    ¬ state.r1.isEmpty = true →
    ops.divmod state.r0 state.r1 = .ok (quotient, remainder) →
    measure (henselEEANextState state quotient remainder) < measure state

/-- Strict well-founded lowering of C++ `polynomial_GCD(F,G,s,t)` over Zp.
It executes the concrete quotient/remainder call at each iteration and then
performs the source `q*s`, subtraction, normalization and six assignments in
order.  The final three inverse-scalings implement source monicization. -/
def __polynomial_GCD_eea_raw_ir (ops : HenselEEARawOps)
    (termination : HenselEEATermination ops) :
    HenselEEAState → RawExec (SparsePolyZp × SparsePolyZp × SparsePolyZp)
  | state =>
      if hdone : state.r1.isEmpty then
        match state.r0[0]? with
        | none => .error .assertionFailure
        | some leading =>
            match ops.inverse leading.2 with
            | .error fault => .error fault
            | .ok inverse =>
                .ok (henselEEAScaleNormalize inverse state.r0,
                  henselEEAScaleNormalize inverse state.s0,
                  henselEEAScaleNormalize inverse state.t0)
      else
        match hrun : ops.divmod state.r0 state.r1 with
        | .error fault => .error fault
        | .ok qr =>
            __polynomial_GCD_eea_raw_ir ops termination
              (henselEEANextState state qr.1 qr.2)
termination_by state => termination.measure state
decreasing_by
  exact termination.decreases state qr.1 qr.2 hdone hrun

def henselEEAInitialState (p : UInt64)
    (left right : SparsePolyZp) : HenselEEAState :=
  let one : SparsePolyZp := #[(UMonomial.mk 0, Zp.ofUInt64 1 p)]
  { r0 := left, r1 := right, s0 := one, s1 := #[], t0 := #[], t1 := one }

end Generated.StrictHensel
