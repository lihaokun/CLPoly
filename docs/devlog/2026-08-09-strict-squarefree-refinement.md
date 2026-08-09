# Strict C++ L1 → Lean L2 SQF refinement

Date: 2026-08-09

## Scope

Close the generated C++ execution of `__squarefree_Zp` against the existing
well-founded L2 `sqfZp`/`yunLoop` model.  The strict path imports the proved raw
public polynomial GCD and supplies its own safe, well-founded definitions for
every generated range-for and recursive loop used by SQF.

## Natural-language proof draft

### Sparse monic helper

1. Read the concrete first array cell, exactly as `f.front()` does.
2. If its stored word is one, return the original array without traversing it.
3. Otherwise evaluate the concrete well-founded `Zp::inv` implementation and
   run the forward mutation loop already proved in the strict GCD module.
4. The loop invariant says that cells before the iterator have been multiplied
   exactly once, later cells are unchanged, monomials never change, and array
   size is preserved.
5. For a canonical nonempty input whose denotation is monic, the first stored
   word is exactly one.  Therefore every monic value produced inside SQF takes
   the real early-return branch; this is proved from the branch comparison,
   rather than replacing the helper with an identity function.

### `p`-th-root extraction loop

1. Traverse the source sparse array from index zero to its actual size.
2. At each iteration read the bounded source cell and append one cell with the
   same coefficient and the source word-level degree division by `p`.
3. Use `source.size - index` as the termination measure.  Since the source is
   immutable and the iterator increases by one, the measure strictly falls.
4. Maintain that the accumulator is its initial prefix followed by the exact
   image of the already visited source prefix.
5. Under canonical degree bounds and the derivative-zero condition, every
   degree is divisible by the characteristic.  Word conversion/division then
   agrees with Nat division, so the final sparse denotation is the L2
   Frobenius contraction used by `sqfZp`.

### Derivative loop

1. Traverse the immutable sparse source in order and skip degree-zero terms.
2. Reproduce `Zp * int64_t` as written by current C++: reduce the nonnegative
   degree into the coefficient modulus, call the generated normalized modular
   multiply, skip a zero product, and append degree minus one otherwise.
3. Use `source.size - index` as the well-founded measure and maintain that the
   output suffix is the exact `filterMap` image of the visited source prefix.
4. The proved generated modular-multiply theorem identifies each emitted word
   with coefficient times degree in `ZMod p`; termwise differentiation and
   additivity then identify the final sparse polynomial with the L2 derivative.
5. Strictly descending source degrees remain strictly descending after
   subtraction by one, while filtering preserves order; generated reduction
   and the explicit zero test establish canonical output.

### Recursive-result exponent loops

1. Traverse the returned factor array using the source range-for index.
2. Append the same sparse factor and the actual `UInt64` product `e * p`.
3. The measure is the immutable source result size minus the iterator.
4. The invariant identifies the destination suffix with the transformed source
   prefix.  A separate no-wrap premise connects the machine product with L2 Nat
   multiplicity multiplication.

### Yun loop and recursive SQF

1. In the derivative-nonzero branch invoke the strict nonempty public GCD
   execution, then the strict sparse exact-division execution for `f/c`.
2. Relate each C++ state `(i,w,c,result)` to the corresponding L2 `yunLoop`
   state: both sparse polynomials are canonical representations of the L2
   values, the result prefix denotes the same factor/multiplicity list, and
   `c` is nonzero.
3. Each iteration invokes the same strict GCD for `y`, exact-divides `w/y` and
   `c/y`, conditionally executes the real monic helper, and appends the same
   factor.  Divisibility and monicity prove both divisions exact.
4. Terminate by the source/L2 measure `degree(w) + degree(c)`; the normalized
   gcd and exact quotient identities give the strict decrease already used by
   L2 `yunLoop`.
5. In derivative-zero and residual-`c` branches extract the concrete `p`-th
   root and recurse on strictly smaller degree, then execute the exponent-copy
   loop.  This matches L2 contraction followed by multiplicity scaling.

## Required evidence

- `lake env lean CLPoly/Refinement/SquarefreeZp.lean`
- targeted Lake build of strict GCD and SQF modules
- source gate for well-founded raw GCD/HGCD
- scan of strict SQF/GCD sources for unchecked placeholders and legacy
  aggregate imports
- `#print axioms` for every exported strict SQF refinement theorem

## 度量

- 终止度量：数组循环使用 `source.size - index`；Yun 使用
  `degree(w) + degree(c)`；递归 SQF 使用输入首项次数加一。
- 当前迭代：10 轮单文件编译—修复。
- 当前状态：公共非空 raw GCD 已闭合；`__upoly_make_monic` 的真实
  分支、`__extract_pth_root` 的前向 append 循环、两个递归结果指数缩放
  循环已经以良基定义形式化，并证明最终数组列表等于对实际已访问前缀
  的精确映射。
- 当前 C++ derivative 的前向循环也已闭合：系数计算调用已认证的生成
  preinverse 模乘，输出保持 canonical，其 `toPoly` 精确等于 L2
  derivative，并且实际 `isEmpty()` 判定与 L2 derivative 为零等价。
