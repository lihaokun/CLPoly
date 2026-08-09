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

### Sparse `pair_vec_div`

1. Preserve the source branch order: reject an empty divisor, handle aliasing,
   clear the destination, return on an empty dividend, then split the
   single-term divisor and general priority-heap paths.
2. In the single-term path traverse every dividend term, run the actual
   univariate monomial divisibility test, and for divisible terms compute the
   coefficient with the generated inverse and normalized modular multiplier
   before appending it.
3. In the general path represent each C++ heap node by its divisor-tail index
   and quotient index.  The heap invariant identifies every live node with the
   product of those two concrete sparse cells; equal monomials form the linked
   bucket consumed by the inner loop.
4. The outer invariant states that `new_v * divisor` plus the unconsumed
   dividend/heap frontier has the same polynomial denotation as the original
   dividend.  Each subtraction bucket preserves this equality, and emitting a
   quotient cell advances the leading frontier.
5. Termination is lexicographic on the finite unconsumed dividend cells, live
   node advances, and remaining output degree.  Exact-division call sites prove
   the final frontier empty and identify the emitted quotient with L2
   `divByMonic`.

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
- `p`-th-root 路径已从数组映射加强为完整语义桥：机器次数转换在
  `int64_t` 次数范围内等于 Nat 除法；derivative 为零推出每个实际
  canonical 项次数可被特征整除；提取结果保持 canonical，且其
  `toPoly` 精确等于 `Polynomial.contract`。
- 递归结果 copy 循环的 multiplicity 语义也已闭合：在明确的
  `UInt64` 无回绕前提下，实际 `e * prime` 的 `toNat` 精确对应 L2
  列表中的 `e * p`；因子多项式保持逐项相同。
- `pair_vec_div` 的 single-divisor source branch 已闭合：逐项 monomial
  divisibility test、generated inverse、generated preinverse multiply 与
  append 均保留；输出保持 canonical，满足 `quotient * divisor =
  dividend`，并在 monic exact-division 前提下精确等于 L2
  `divByMonic`。一般 divisor 的 priority-heap branch 仍待闭合。
- 一般 divisor 路径已经开始按源布局建模：每个 non-leading divisor
  cell 对应一个 `VHC` 节点，quotient/divisor 指针转成有界索引，尚未
  激活的未初始化 `mono` 明确表示为 `none`；节点分配循环及其精确
  长度/列表定理已经闭合。`reset_h` 的单节点激活也已实现为 checked
  raw→safe 步骤，并证明只写目标节点、保持节点数组长度。下一检查点
  中，`VHC_insert`/`VHC_extract` 已进一步形成可编译的 checked 执行：
  插入保留空 heap、root 等次桶、新最大 root、内部 anchor 等次桶与
  普通上浮五条分支；extract 保留末槽 sentinel、左右 child 比较、
  下沉复制和最终逻辑 pop。父节点搜索、两种上浮与下沉均有严格下降
  的 Nat 度量。heap pointer、max-order、等次 `next` 链及节点乘积含义
  的组合不变式已经显式定义；下一步证明 insert/extract 保持这些
  不变式，再接 inner bucket consumption。
- inner `while (heap[0] != nullptr)` 已开始闭合：单节点执行使用 generated
  `nmod_mul` 后接 generated `nmod_sub` 实现真实 `submul`，其 ZMod 语义
  已证明为 `k - qCoeff * divisorCoeff`；节点随后严格按源码二分为推进
  `quotientIndex` 并压入 `lin`，或耗尽并增加 `reset_h`，且证明两种
  进展恰有一种发生。完整 `next` 链使用未消费节点有限所有权集作
  良基度量，每步删除当前节点并检查环，不使用 fuel。下一步先给
  sift/extract 建立保持长度/严格缩短定理，再闭合相同 monomial 的外层
  bucket 循环和 `lin` 逆序重插。
- 上述 heap 控制流现已闭合：手工归纳证明 sift 仅改写槽位、保持数组
  长度，进而证明每次成功 extract 恰好把逻辑 heap 大小减一，并形成
  proof-carrying checked extract。相同 monomial 的外层 bucket 循环因此
  直接以实际 `heap.size` 良基递归。迭代末尾的 `lin[--lin_size]` 逆序
  重插和新 quotient cell 触发的 `--reset_h` 节点激活/插入也已按源码
  顺序实现，分别以 `lin.size` 与 `reset_h` 终止。下一步组合一次完整
  outer iteration，并建立 quotient-product frontier 不变式。
- 一次完整 outer iteration 现已组合为 checked 执行：严格按 source
  tie-breaking 在 dividend cell 与 heap root 中选择 frontier（相同次数
  先推进 dividend iterator），消费等次 bucket，执行非零与 monomial
  divisibility 分支，以 generated inverse/multiply 生成 quotient cell，
  仅在真实 `value != 0` 分支 append，随后依次运行 `reset_h` 激活和
  `lin` 逆序重插。已证明 frontier 选择只会保持或恰好推进一次实际
  dividend iterator。下一步建立本次迭代后最大 frontier 次数严格下降
  的不变式，用它定义完整 outer while 的良基递归。
- 审计源码后修正了 monic helper 的逆元执行：不再调用旧对象模型的
  inverse，而是构造 `Zp` 结果于已经认证的 generated `inv_prime` word
  execution；monic 早退分支仍按真实 stored-word comparison 执行。
