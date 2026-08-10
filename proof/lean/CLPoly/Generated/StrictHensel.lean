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

end Generated.StrictHensel
