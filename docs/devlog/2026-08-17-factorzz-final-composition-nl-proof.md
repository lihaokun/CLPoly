# FactorZZ strict final composition: natural-language proof draft

## Goal

Construct one concrete `Generated.StrictFactorZZ.FactorZZRawOps` whose four
fields execute the already-generated strict C++ L1 entries for prime
selection, Hensel lifting, van-Hoeij recombination, and Zassenhaus
recombination.  Prove that a successful execution of the generated
`__factor_squarefree_primitive_ZZ_raw_ir` returns the same physical array that
satisfies `FactorZZCorrect`.

No field may return a proposition, a chosen factorization, an L2 result, or a
semantic fallback.  Runtime readiness assumptions may provide arithmetic
storage/workspaces and universally quantified safety invariants, but may not
prescribe any intermediate or final output.

## Concrete operation bundle

1. `selectPrime` closes the actual generated `__select_prime_raw_ir` over the
   concrete GMP prime iterator, modular reduction, derivative/GCD, DDF, EDF,
   fixed initial RNG state, and the per-prime physical workspace provider.
2. `henselLift` checks the runtime word for primality, obtains only the dense
   arithmetic/workspace object indexed by that word, and executes the actual
   generated `__hensel_lift_upoly_raw_ir` with the concrete division
   termination certificate.  A non-prime word returns the source-shaped
   assertion fault; it cannot produce polynomial data.
3. `vanHoeijRecombine` checks the source degree precondition and executes
   `Generated.StrictRecombine.__vanhoeij_recombine_raw_ir` with
   `concreteVanHoeijRawOps` and `concreteVanHoeijTermination`.
4. `zassenhausRecombine` is literally
   `Generated.StrictRecombine.zassenhausRecombine
   concreteZassenhausTermination`.

## Proof of the top-level controller

Let the successful prime-selection output be `selection`.

1. Expand only the generated `__factor_squarefree_primitive_ZZ_raw_ir` entry.
   The source input hypotheses rule out its assertion branch.
2. Apply the existing strict select-prime refinement theorem to the exact raw
   call.  This yields `SelectionCorrect source selection`, hence primality,
   leading-coefficient survival, squarefreeness of the modular image, a
   product association with `selection.factors`, and irreducibility/monicity
   of every stored modular factor.
3. If the literal source guard reports `selection.irreducible` or at most one
   modular factor, the generated output is exactly `#[f]`.  Use
   `SelectionCorrect`'s one-factor/irreducible consequence plus primitivity of
   the integer input to prove `Irreducible (toPoly f)`; the singleton product
   is definitionally associated with the source.
4. Otherwise the source has at least two modular factors and executes the
   generated `__lll_factorize_raw_ir`.  Expand its successful heuristic call,
   first Hensel call, and first van-Hoeij call.  All values are the unique
   physical outputs of those calls.
5. Instantiate the Hensel entry theorem using the workspace indexed by
   `selection.prime` and the universally quantified entry invariant.  Derive
   the exact Hensel output certificate, including equal cardinality with the
   selected factors, canonical heads, prime/prime-power product association,
   pairwise coprimality, and full-precision bound when target zero is used.
6. Apply the generated van-Hoeij entry theorems to obtain:
   (a) output size at most lifted size, and
   (b) association of the integer source with the product of the exact
   physical output array.
7. Classify the literal controller branches.
   - If a low-precision output is accepted, the signed-size comparison plus
     the size upper bound proves equality with the lifted cardinality.
   - If low precision is rejected, execute the actual second Hensel and
     van-Hoeij calls.  At full precision, a smaller result triggers the
     literal Zassenhaus call; otherwise the size upper bound makes the result
     equal-cardinality.
   - If the first call was already full precision, the same safety predicate
     either executes Zassenhaus for a smaller result or leaves an
     equal-cardinality result.
8. In every actual Zassenhaus branch, apply
   `selectionHensel_zassenhausRecombine_refines_FactorZZCorrect` to that exact
   call and output.
9. In every equal-cardinality van-Hoeij branch, prove physical member quality:
   every accepted validation factor and remaining quotient has nonzero
   leading coefficient and is a nonunit after reduction.  Combine these facts
   with the Hensel modular irreducibles, equal cardinality, and the exact
   product association using
   `factorArrayIrreducible_of_hensel_cardinality`.  This proves irreducibility
   of the same array returned by van-Hoeij, not of an existentially chosen
   replacement.
10. Pair the product association from step 6 with member irreducibility from
    step 9 to obtain `FactorZZCorrect`.

## Required missing lemmas

The existing sources already provide every item above except the following
entry-level packaging for the equal-cardinality branch:

1. `__vanhoeij_recombine_raw_ir_members_nonunit_mod` for the exact returned
   array, derived through `vanHoeijLoop` from candidate validation and the
   concrete Zassenhaus fallback.
2. `__vanhoeij_recombine_raw_ir_leading_mod_ne_zero` for the exact returned
   array, derived through the same execution trace.
3. A composition lemma applying the preceding facts,
   `__vanhoeij_recombine_raw_ir_product_associated`, the output-size theorem,
   and `factorArrayIrreducible_of_hensel_cardinality`.

These are execution invariants.  They may quantify over the raw output of a
successful call, but must not contain an existential factorization witness or
invoke UFD factorization to replace the generated result.

## Progress after the draft

The first representation bridge needed by the loop invariant is now proved:

- `gatherActive_members_source` follows the literal checked gather loop and
  proves that every returned physical polynomial is a member of the original
  lifted array;
- `gatherActive_irreducible` transfers modular irreducibility from that
  original Hensel array to the exact gathered subarray used by candidate
  validation and the internal Zassenhaus fallback.

Both lemmas compile in `CLPoly.Refinement.Recombine` and are included in the
refinement axiom audit.  They do not assume valid indices in advance: a
successful raw execution itself proves that both source guards were passed.

## Audit of the draft

- Mathematical correctness: equal cardinality plus association to the same
  list of modular irreducibles rules out grouping; the full-precision smaller
  branch is handled by exhaustive Zassenhaus.
- No skipped controller branch: assertion, irreducible/singleton, low/full
  precision, equal/smaller cardinality, and both Hensel calls are enumerated.
- Lean path: every named dependency exists; the three missing entry invariants
  are isolated before the top theorem.
- Termination: all executed loops use existing structural, natural,
  determinant-potential, active-set, or combination-rank well-founded
  certificates.  No new recursion is required by the composition theorem.
- Boundary cases: the theorem assumes canonical, nonempty, primitive,
  degree-at-least-two input; the generated assertion branch is therefore
  unreachable.  The selected one-factor branch and both possible
  full-precision paths are explicit.

## 度量

- 耗时：约 0.8 小时（现状审计、依赖定位和证明草稿）。
- 迭代：尚未进入 Lean 编译迭代。
- Lean 新增/修改行数：约 65 行（active-gather provenance and transfer）。
- 对应 C++ 行数：顶层控制器约 80 行。
- 放弃的方案：将 van-Hoeij 结果映射到某个抽象不可约分解；该方案不验证返回的物理数组，明确禁止。
