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
- 完整 outer while 与统一 dispatch 已形成可编译执行。while 每轮读取
  真实 frontier，并用“上一实际 frontier 的严格次数上界”作为良基
  参数；递归参数更新为本轮真实次数，不是迭代计数或 fuel。若次数不
  下降则返回 checked assertion failure，后续语义证明必须从 canonical
  frontier/heap 不变式排除该分支。一般分支入口以 dividend 首项次数
  加一建立初始上界，并创建真实 `divisor.size - 1` 节点/reset_h。
  dispatch 现按源码覆盖 zero divisor、zero dividend、single divisor 和
  compression-false general heap；下一步证明 decrease guard 恒成功及
  quotient-product frontier 守恒，尚未据此宣告一般除法精化完成。
- decrease guard 的第一层语义现已闭合：显式 frontier-bound 不变式
  同时约束所有未消费 dividend cell 与所有活跃 heap node 的次数；已
  证明真实 selector 的结果严格低于该上界。canonical sparse dividend
  的严格降序链进一步证明每个实际数组 cell 的次数不超过首项，因此
  空 heap 的一般分支初态确实满足 `head.degree + 1` 上界，而非把初始
  guard 当作假设。剩余工作是证明 outer iteration 把此不变式收紧并
  保持到本轮真实 frontier 次数，同时加入 quotient-product 守恒部分。
- iteration preservation 的 dividend 半边现已闭合：从 canonical chain
  得到任意两个实际数组索引 `i < j` 时后项次数严格更小；据此证明
  selector 无论选择 dividend（iterator 前进一步）还是选择更大的 heap
  root，返回 iterator 后的所有 dividend cell 都严格低于本轮 frontier。
  heap 半边采用更强的全节点块不变式：所有已激活节点无论暂居 heap、
  `next`、`lin` 或 retired suffix，其 monomial 都低于 frontier；该性质
  立即推出 heap-side bound，并已证明单独修改 `next` 指针保持它。下一
  步证明 insert/extract 保持全节点 monomial，再处理 consume/activate
  唯一会重算 monomial 的两处。
- 已进一步把 selector 的 dividend 结论封装为完整保持定理：返回的
  `frontier.dividendIndex` 之后每个实际 cell 都严格低于返回 frontier
  次数，覆盖 dividend-selected、heap-selected 与 dividend-exhausted
  三类真实分支。全节点 bound 的 `setNext` 保持证明也已闭合到具体
  `Array.set` 的同址/异址两种读取。`VHC_insert` 的组合保持仍待证明；
  其执行定义已从 do-notation 等价重写为逐调用的显式 raw match，五条
  源码成功/错误路径均直接可见，后续保持证明不再依赖 monad 大展开。
- `VHC_insert` 的五条成功分支现已全部证明保持全节点 frontier 上界，
  并组合成不要求调用者预先分类分支的统一保持定理；`lin` 的逆序重插
  循环也已用其实际数组长度作良基归纳，证明任意多次真实 insert 均保持
  该上界。下一步只剩两类会重算节点 monomial 的状态变化：消费节点时
  推进 quotient 指针，以及 append 新 quotient 后的 `reset_h` 激活。
- 消费节点推进 quotient 指针的核心语义现已闭合：由 quotient canonical
  严格降次与节点的真实 `quotientIndex`/`divisorIndex` 表示，证明重算后的
  乘积 monomial 严格小于原节点；据此 `pairVecDivVHCConsumeNode` 保持全
  节点 frontier 上界，同时保持每个活跃节点确实表示对应 quotient/divisor
  数组单元乘积。接下来把这两个性质沿良基 `next` chain 抬升，再处理
  `reset_h` 激活。
- 在沿 chain 抬升时重新逐行核对 `basic.hh`，发现并修正了 exhausted
  分支的真实性缺口：源码的 `++v1_ptr` 即使到达 `v_size` 也会留下
  one-past-end 指针，旧安全模型却保留了旧 quotient 索引。现在该分支
  保存递增后的索引，并把直到 `reset_h` 激活前不可观察的陈旧
  `mono`/`next` 标成 `none`；这既匹配源码可观察行为，也使无效指针状态
  无法被误读。可达 `next` 链已建立独立的 Finset 所有权不变量，允许
  真实未初始化 reset 后缀，并已证明单节点更新保持其剩余链；完整 chain
  执行现保持 frontier 上界和活跃节点 denotation。
- `reset_h` 单节点激活的两项局部保持也已闭合：真实 activate 写入的新
  monomial 在给定本轮界时保持全节点 frontier 上界，并重新建立该节点
  对 quotient/divisor 单元乘积的 denotation。另已从 divisor canonical
  严格尾项降次证明关键算术桥：新 quotient 次数为
  `frontier - divisorLead` 时，它与任意非首 divisor 项的乘积次数严格低于
  frontier。下一步把 reset 前缀的 one-past-end 指针关系沿激活循环保持。
- reset 前缀关系现已显式化为 `PairVecDivVHCResetReady`：前缀中每个节点
  的 quotient 索引等于旧 quotient 长度、divisor 索引保持 `i + 1`，且
  product monomial 尚不可观察。节点初始化数组已证明满足该关系（包括
  每个具体数组索引的精确初值），为后续证明 exhausted 扩张与 activate
  收缩同一个真实 `reset_h` 前缀提供了基例。
- 源码中永久不变的 `v2_ptr` 现对应独立不变量
  `PairVecDivVHCNodeDivisorIndicesFixed`，精确陈述节点槽 `i` 始终指向
  divisor 尾项 `i + 1`。该性质已从节点分配证明，并贯穿 consume 的推进/
  exhausted 两支、`setNext`、五分支 insert、`lin` 良基重插、activate，
  以及完整 `ActivateReset` 良基循环。另已证明在 consumed node 恰为当前
  `resetH` 时，真实 exhausted 分支把 reset-ready 前缀从 `resetH` 扩为
  `resetH + 1`；剩余关键点是从 heap/bucket 次序证明这一节点编号条件。
- exhausted 对源码裸 `++reset_h` 的行序依赖现已提升为 checked raw guard
  `nodeIndex = resetH`。这不是额外执行预算：合法 C++ 状态仍须由 heap
  枚举不变量证明 guard 恒成立；但任何成功 safe 执行从此无法绕过该
  义务。利用 guard、活跃节点不可能位于 reset 前缀以及永久 `v2_ptr`
  身份，已证明任意成功 consume node 保持 `ResetReady`，并进一步沿有限
  `next` chain 良基递归抬升；完整 chain 结果现在同时保持 frontier、
  denotation、节点身份和 reset 前缀四项不变量。
- reset 前缀的收缩侧也已闭合到完整源码循环：建立了“缩短前缀并在前缀
  外写一个节点”的通用保持引理，证明 activate 写入 `resetH - 1` 后剩余
  前缀仍 ready，随后真实 insert 的 `next` 写入也不影响该前缀；以实际
  `resetH` 良基递归得到 `pairVecDivVHCActivateReset` 成功后
  `ResetReady 0`。下一步把单节点已有的 frontier/denotation 保持与此
  收缩证明同步组合，再抬升 equal-degree heap bucket。
- `ActivateReset` 现已从单纯前缀终止加强为四项组合保持定理：在新
  quotient 尾项次数等于 `frontier - divisorLead` 的真实 append 关系下，
  每次 activate 利用 canonical divisor 尾项严格降次证明新 heap product
  低于 frontier；随后的 insert 同时保持 frontier、denotation、永久
  divisor 身份和剩余 reset 前缀。完整良基循环最终给出四项性质及
  `ResetReady 0`，不再把语义保持留给调用端假设。
- `next` chain 的四项组合定理现已封装到真实
  `pairVecDivVHCConsumeRootBucket` 调用边界：从非空 heap 的实际 root 和
  `Finset.range nodes.size` 所有权集合执行，不再要求上层展开 chain
  实现。下一步需要加强 heap ownership，使 `VHC_extract` 后的新 root
  自动取得未被当前 bucket 修改的合法 chain，而不是由调用端逐轮提供。
- bucket safe 结果现增加纯证明可见的剩余 `unvisited` 集合，不改变任何
  C++ 可观察字段。已沿良基 consume chain 证明：结果集合包含于入口集合，
  且其中每个节点在结果数组中逐点等于入口数组；root-bucket 边界也已
  得到对应定理。另证明 `ChainValid` 对所有权集合内逐点相等封闭。这样
  heap ownership 下一步只需证明其他 bucket 的不交集合仍包含在剩余
  `unvisited` 中，即可把其 chain 合法性传过当前 bucket 的节点修改。
- chain ownership 已从“合法可达”加强为精确 `PairVecDivVHCChainOwns`：
  沿 `next` 走完时 owner 集合必须恰好耗尽，并证明其所有索引都落在真实
  节点数组范围内。核心隔离定理现已闭合：任意与当前 root owner 不交的
  另一 chain，在真实 root-bucket 消费后，其 owner 全部仍属于结果
  `unvisited`，对应节点逐点未变，因而 `ChainOwns` 自动迁移到结果节点
  数组。下一步把每个 heap slot 配上两两不交 owner，再证明 extract 只
  重排这些 slot/head，不改变 owner 对应关系。
- `VHC_extract` 的真实 sift-down 执行现已补齐两个关键语义引理：
  结果每个活动 head 都有旧 heap slot 来源；若旧 root head 唯一，
  extract 后它在整个活动 heap 中完全消失。证明逐步跟踪了 C++
  sift-down 的 child-copy 和 saved-last-node 写入，没有把 extract 替换成
  抽象排序。剩余局部缺口是证明这些写入不复制存活 head，以得到
  slot 唯一性的 extract 保持。
- 全 heap 的 bucket 所有权现已用 `PairVecDivVHCHeapChainOwnership`
  显式建模：每个 head 精确拥有一条会耗尽 owner 集合的 `next`
  chain，head slot 唯一，不同 head 的 owner 两两不交。由此已将
  root-bucket 的局部隔离定理提升到全 heap：真实消费 root chain
  后，所有其他活动 head 仍精确拥有原 chain。下一步仅需把
  extract 的 slot 唯一性保持与已有 provenance/root-removal 引理组合。
- extract 的 slot 唯一性现已闭合。证明使用真实 sift-down 路径
  上的“前缀唯一，仅当前 hole 可临时重复”不变式：child 被复制
  到 hole 后将例外转移到 child，最后 saved-last-node 写入 hole，
  `pop` 删除 sentinel 位置后恢复全前缀唯一。该定理已与 head
  provenance、root 彻底移除、其他 chain 不变组合，得到单次
  root-consume-plus-extract 对完整 heap chain ownership 的保持。下一步沿
  heap-size 严格下降将它提升到整个 equal-degree 良基循环。
- 完整 `pairVecDivVHCConsumeEqualDegree` 的 heap chain ownership 保持已
  闭合。证明直接对真实活动 heap size 做强归纳，每次递归使用
  `pairVecDivVHCExtractChecked` 包装的原始 extract 执行以及严格
  `size` 下降证据。因此该循环现保持每个存活 bucket 的精确
  owner、head 唯一性和 owner 两两不交；未引入 fuel 或替代循环。
- equal-degree 良基循环现亦同步保持四类后续除法所需不变式：
  所有活动节点的 frontier 次数界、节点对具体 quotient/divisor
  product cell 的 denotation、节点到 divisor 行的固定对应，以及
  `reset_h` 前缀的一步越界状态。root chain 合法性不再作为额外
  假设，而是由全 heap 精确 ownership 及 owner 落在节点数组范围自动
  推出。
- 为处理 activation/reinsert 会改变 bucket 链拓扑的事实，全局
  不变式已提升为“存在当前 owner 映射”的
  `PairVecDivVHCHeapChainsOwned`；空 heap 初始化和整个 equal-degree
  循环对它的保持已闭合。同时已证明真实 `VHC::next` 写入的
  ownership 效果：将一个 fresh node 链到旧 chain 前方，新 head 精确
  拥有 `insert nodeIndex tailOwner`。这是下一步证明真实
  `VHC_insert` 保持 existential heap ownership 的核心链接引理。
- 两条真实指针复制路径 `pairVecDivVHCBubble` 和
  `pairVecDivVHCBubbleBelow` 现已证明保持 head provenance、数组长度与
  唯一性。唯一性证明跟踪“只有当前 hole 可临时与 parent
  重复”和“fresh head 只在当前 hole”两个前缀不变式，直到最终
  head 写入。另补充了一般 heap-reordering ownership 迁移定理以及
  `SetNext` 对不交 chain 的保持。剩余 insert 工作是在等次合链
  与 fresh singleton 入堆分支中构造新 owner 映射。
- fresh singleton 入堆的 owner 映射现已构造完成：新 head 的 owner
  为 `{newNode}`，旧 head 保留原 owner。证明同时要求并使用了新
  node 不在旧 heap 和任一旧 owner 中的 freshness，因此新 singleton
  与旧 chain 两两不交。该 ownership 已沿真实 `Bubble` 和
  `BubbleBelow` 执行迁移到最终 heap，两条新 bucket 入堆分支的
  表示不变式核心已闭合。剩余 `VHC_insert` 局部缺口是 root/anchor
  等 monomial 时的 owner-union 合链。
- root/anchor 等 monomial 的 owner-union 合链现已闭合：旧 head 被唯一
  fresh head 替换，新 owner 为 `{newNode} ∪ oldOwner`，其他 bucket
  owner 不变；证明包括新链精确性、替换后 head 唯一性与所有
  跨 bucket 不交性。在此基础上，完整 `pairVecDivVHCInsert`
  的所有真实分支已统一证明 existential heap ownership 保持，前提
  是 activation/reinsert 提供的新节点 freshness。下一步从旧节点
  `mono = none` 和现有 exact-chain ownership 推出该 freshness，并将定理接入
  activation/reinsert 良基循环。
- freshness 现已从表示不变式自动推出：exact chain 的每个 owner
  成员都必然是 `mono = some _` 的活动节点，因此
  `reset_h` 前缀中 `mono = none` 的节点不可能已在 heap 或任一
  活动 owner 中。利用该结论，完整 `pairVecDivVHCActivateReset`
  循环已证明 existential heap ownership 保持：每步执行真实
  activate 写入，再执行真实 `VHC_insert`，并以 `resetH - 1` 作为
  良基下降。剩余表示层缺口是为 `lin` 反向重插证明其节点
  两两唯一且不在当前 heap owners 中。
- 为 `lin` 反向重插建立了
  `PairVecDivVHCHeapChainsOwnedAway protectedSet`：除精确 heap
  ownership 外，显式记录 protected 节点不是 heap head，且与每个
  bucket owner 不交。已证明空 heap 初始化、protected 集合缩小及
  从 protected 中取出一个活动节点时的 `VHC_insert` freshness/
  ownership 结论。下一步将 insert 后的 away 集合精确更新为
  `protectedSet.erase newNode`，以便对 `lin.pop` 做良基归纳。
- 审计源码后修正了 monic helper 的逆元执行：不再调用旧对象模型的
  inverse，而是构造 `Zp` 结果于已经认证的 generated `inv_prime` word
  execution；monic 早退分支仍按真实 stored-word comparison 执行。
- `lin` 反向重插的表示不变式现已严格闭合。完整真实 `VHC_insert`
  的空堆、新 root 上浮、anchor 下方上浮、root/anchor 等次数合链分支，
  均会把待处理 protected 集合精确缩小到尚未重插的节点，并保持它们
  与所有 heap chain owner 不交。新增 `PairVecDivVHCLinReady` 记录
  `lin` 节点无重复且全部处于活动状态；证明每次 `lin.pop` 删除的恰是
  当前末节点，其他节点的数组内容不受 `SetNext` 写入影响。最终对真实
  `pairVecDivVHCReinsertLin` 按 `lin.size` 作良基递归，得到重插结束后的
  完整 existential heap-chain ownership；未使用执行步数参数、规格结果
  或替代算法。下一步把 equal-degree 消费产生的 `lin`/heap 状态直接
  证明满足 `Away + LinReady`，从而接通 outer iteration。
- equal-degree 输出不变量的第一层现已闭合到真实单节点消费：若当前
  chain 节点尚未进入 `lin`，`pairVecDivVHCConsumeNode` 的 advance 分支
  对真实 `nodes.set` 与 `lin.push` 同时保持 `LinReady`，exhausted 分支
  则证明 `lin` 中既有活动节点不受该写入影响；两分支还给出输出
  `lin` 集合只可能增加当前节点的精确上界。下一步沿 exact chain owner
  的良基消费归纳提升到整个 root bucket。
- `Away + LinReady` 已从单节点提升到完整 equal-degree 良基循环。
  `pairVecDivVHCConsumeChain` 证明输出 `lin` 只来自入口 `lin` 与当前
  exact chain owner；因此节点不重复且保持活动。root bucket 消费后，
  真实 `VHC_extract` 删除当前 head，并借助不同 bucket owner 两两不交，
  证明存活 heap chain 与扩展后的 `lin` 仍完全隔离。最终按真实 heap
  size 严格下降提升到 `pairVecDivVHCConsumeEqualDegree`，其结果直接满足
  后续 reverse reinsertion 所需的 `HeapChainsOwnedAway + LinReady`。
  下一步将此组合定理接入 outer iteration 的 activate/reinsert 分支。
- outer iteration 的 activation 隔离缺口已闭合。`VHC_insert` 的完整
  away 保持核心被推广为接受由表示不变量证明的显式 freshness；
  protected-node 重插仍由原包装定理自动推出 freshness。对真实
  `reset_h` activation 循环，`mono = none` 证明待激活节点不属于活动
  `lin`，exact ownership 证明它不在 heap/owner 中；随后真实 activate
  写入和 `VHC_insert` 保持同一个 `lin` protected set。按 `resetH`
  良基递归，最终同时保留 `HeapChainsOwnedAway` 与 `LinReady`，可直接
  交给已经闭合的 reverse reinsertion 定理。
- 完整 outer iteration 的 heap-chain ownership 现已闭合。商项
  emission 被提取为与原始 C++ 分支一致的 `pairVecDivVHCEmit`：
  它显式携带 divisor 非空证明，执行真实 modular inverse/multiply、
  quotient push 与 `reset_h` activation，不用规格结果替换计算。
  `pairVecDivVHCOuterIteration_preserves_heapChainsOwned` 现在直接组合
  frontier selection、equal-degree consumption、emit/activation 与良基 reverse
  `lin` reinsertion，得到一次完整外层循环体保持 existential exact-chain
  ownership。下一步将这个保持定理提升到完整良基 outer loop，
  然后闭合 general division 语义。
- outer iteration 的其余表示不变式已与 ownership 同步接通。
  emit/activation 证明保持活动节点次数上界、节点对新 quotient
  的 denotation、固定 divisor index 以及 `reset_h` ready 前缀；其中
  显式证明了 quotient `push` 不改变所有旧项的可寻址性。reverse
  `lin` reinsertion 也已按 `lin.size` 良基递归保持这四类不变式；
  对于非零 `reset_h`，每个待重插节点在 ready 前缀之外是由
  `lin` 活动性与 ready 前缀的 inactive 性质直接推出，不是附加假设。
  完整单次 outer iteration 现返回 ownership 与全部节点不变式。
  剩余 loop 缺口是将次数上界从旧 `degreeLimit` 收紧到实际
  selected frontier，以供下一次良基调用。
- 完整 loop 所需的 global frontier-bound 骨架已建立。新不变式
  `PairVecDivVHCRemainingDividendBelow` 记录当前 source index 之后的
  所有 dividend 项都低于固定全局上界；selector 只会保持或递增
  source index，因此该不变式可直接传递。对 heap 侧，exact-chain
  ownership 证明每个可见 head 是真实活动节点，再由
  `AllActiveNodesBelow` 得到 head 上界；两者组合生成 selector 所需
  `FrontierBelow`，无需假定所有活动节点都已在 heap 中。单次
  outer iteration 现已同时传递 remaining-dividend bound、ownership 和
  全部节点不变式；canonical 非空 dividend 的初始上界也已闭合。
- quotient emission 所需的 coefficient reducedness 已沿真实 accumulator
  执行闭合。selector 返回 dividend 存储系数或零，其范围分别由
  dividend canonicality 与模数正性推出。每个消费节点执行真实
  generated modular multiply 与 subtraction；multiply 的 preinverse 正确性给出
  product 范围，source subtraction 的 range theorem 给出新 accumulator
  `< p`。该性质已按 owner-set card 良基递归提升到完整 next-chain，
  再按 heap size 良基递归提升到完整 equal-degree 循环。下一步用它
  证明 nonzero emitted quotient coefficient 是 canonical 系数。
- quotient emission 的 canonicality 已闭合到真实 emit 执行。generated
  modular multiply 在已证 accumulator `< p` 的前提下给出 emitted
  value `< p`；`value ≠ 0` 是 source 真实 push 分支的 guard；已有
  quotient 项次数严格更高时，`canonical_push_lower` 证明真实
  `Array.push` 保持 canonical sparse representation。其他 emit 分支不改变
  quotient。同时已证明：对任意严格低于当前 selected frontier 的
  下一次数，emit 后的所有 quotient 项仍严格位于其上。剩余工作是
  将该次数与下一次真实 selector 执行关联。
- canonical quotient 已接入完整单次 outer iteration。真实 selector
  从 canonical dividend 得到初始 coefficient 范围，equal-degree 消费循环
  保持 accumulator reducedness，emit 再由 preinverse-configured generated multiply
  与 quotient-degree invariant 证明输出 quotient canonical。reverse reinsertion
  不改变 quotient。此外新增完整 iteration 级别的 look-ahead 定理：
  给定下一 selector 次数严格更低，重放真实 consume/emit/reinsert 执行
  即可证明输出 quotient 的所有项仍严格高于下一 frontier。这是
  完整良基 outer loop 递归传递 canonicality 的最后局部桥。
- 完整良基 outer loop 的 quotient canonicality 已闭合。成功的非终止
  loop 执行现可分解为：真实 selector guard 严格下降、一次真实
  outer iteration 成功、以及真实递归 loop 成功。另一执行桥证明
  任何成功 loop 中实际 selector 的次数必严格小于该调用的
  `degreeLimit`。主定理按 `degreeLimit` 做强归纳：用当前 iteration
  证明传递 canonicality/ownership/节点不变式，再从下一次真实
  selector guard 推出 quotient-degree invariant。因此递归前提由成功执行
  本身建立，未引入 fuel、oracle 或“假定下一轮正确”。
- general division 的第一个系数语义原子已闭合。从一次成功的
  `pairVecDivVHCConsumeNode` 执行中，现可推出 node、quotient 和
  divisor 三个真实访问均在界内。在 prime/preinverse 配置与 canonical
  quotient 前提下，generated `submul` 正确性给出精确等式：
  输出 accumulator 的 `ZMod p` 值等于输入值减去该 node 当前
  `quotientIndex` 与固定 `divisorIndex` 指向的两个存储系数之积。
  该定理直接描述真实节点字段和真实执行，不是事后把差值定义为
  “语义”。下一步将其沿 next-chain 消费轨迹求和。
- next-chain 的 coefficient 求和语义已闭合。新建立的
  `PairVecDivVHCConsumeTrace` 与真实 `ConsumeChain` 一步一一对应：
  done 构造子对应 null head；step 构造子携带真实 `ConsumeNode`
  成功等式、owner-set membership、三个数组边界证明，并记录节点在
  更新前实际读取的 quotient/divisor coefficient pair。按 `unvisited.card`
  良基递归证明每次成功 chain 执行都产生该轨迹；再按轨迹归纳得到
  `result.coefficient = initial coefficient - Σ(qᵢ·dᵢ)` 的精确 `ZMod`
  等式。下一缺口是显式保持 next-edge 的 equal-monomial 性质，以证明
  该求和中每个乘积都属于当前 frontier 次数。
- next-chain 中每个真实乘积的 frontier 次数归属已闭合。新增的
  `PairVecDivVHCChainAtDegree` 按实际 `next` 指针和 owner 集合良基递归，
  并由 `PairVecDivVHCChainValid + PairVecDivVHCNextValid` 以及真实 head
  单项式推出整条链同次数。单节点消费只改写当前节点，因而保持实际
  `next` 指向的尾链次数；结合 node denotation 的保持定理，对真实
  `PairVecDivVHCConsumeTrace` 归纳证明每个记录的 coefficient pair 都确实
  来自 quotient/divisor 中可寻址的存储项，且两项次数之和等于当前
  frontier。这里没有构造规格乘积列表，也没有用输出倒推输入。
  下一缺口是把该定理提升到 root bucket/equal-degree 循环，并补齐执行
  过程中 `PairVecDivVHCNextValid` 的保持，使每个被取出的 bucket 都能建立
  同次数前提。
- root bucket 的真实 coefficient 语义现已组合闭合。给定 heap root 的
  实际 `pairVecDivVHCMono` 读取、exact chain owner、入口 `NextValid` 与
  node denotation，先把 owner 链扩张到真实 node-array range，再证明整条
  root 链与 selector 次数相同；随后重放 `ConsumeRootBucket`，同时得到
  实际消费 trace、`ZMod` 累加器求和等式，以及 trace 中每个存储乘积的
  frontier 次数归属。审查中确认全局 `NextValid` 在消费中并不应被保持：
  C++ 会让已消费/待重插节点暂时保留陈旧 `next`。后续提升将只给仍由
  heap 拥有且彼此不交的链携带同次数性质，忠实覆盖这种中间态。
- heap-owned bucket 的同次数不变量已贯穿完整 equal-degree 循环。
  `PairVecDivVHCHeapChainsHomogeneous` 只量化真实 heap slot 及其 exact
  owner chain；入口可由全局 `NextValid` 一次性建立，但消费之后不再要求
  已移入 `lin` 的暂态节点满足陈旧 `next` 关系。root chain 消费利用 owner
  两两不交证明其他 bucket 的每个节点读取完全未变，真实 `VHC_extract`
  又证明新 heap 的每个 head 都来自旧的非 root slot，因此 ownership 与
  chain degree 同时保持。最后按 checked extract 给出的实际 heap-size
  严格下降完成 `ConsumeEqualDegree` 归纳。下一步用该不变量在每次循环
  迭代调用 root bucket 求和定理，并把各 bucket 的实际乘积 trace 串接。
- 完整 equal-degree accumulator 语义已闭合。新增 denotation 的独立
  chain/root 保持定理，使递归调用无需捎带无关表示假设；owner-chain
  次数性质可单调扩张到真实 node-array range。主定理按 checked extract
  的实际 heap size 强归纳：每轮从 `ConsumeRootBucket` 的真实执行 trace
  取得本 bucket 的存储乘积列表，从递归执行取得后续列表，按执行顺序
  用 `List.append` 拼接，并证明最终 coefficient 的 `ZMod` 值精确等于
  初值减去拼接列表的乘积和。列表中的每个项均保持“来自真实
  quotient/divisor 存储项且次数和等于当前 degree”。空 heap 与新 root
  次数不同的真实退出分支都产生空列表。下一步把该 coefficient-level
  等式与 sparse polynomial coefficient/evaluation 表示连接，形成 general
  division 的多项式等式。
- dividend 指针的次数分割不变量已建立。新的 execution-indexed
  `PairVecDivVHCFrontierSource` 证明 selector 成功时只能来自两种真实
  C++ 读取：读取当前 `dividend[dividendIndex]` 并将指针加一，或读取
  实际 heap root monomial 并以零初始化 coefficient、保持指针不变；
  不允许凭规格任意选择次数。`PairVecDivVHCConsumedDividendAbove` 与已有
  `RemainingDividendBelow` 配对，记录指针前的已处理项均不低于下一
  frontier、指针后的未处理项均低于旧严格界。selector 的真实 decrease
  guard 保持此前缀性质，且 outer iteration 的实际返回 index 被证明精确
  等于 selector 返回 index，因此该不变量已经接通一次完整循环体。
  下一步利用 canonical 严格次数顺序，把这一前缀/后缀分割转成
  `toPoly.coeff frontier.degree` 的精确读取定理，再与 equal-degree 乘积和
  合并。
- selector coefficient 已接到真实 dividend 多项式系数。heap-source
  构造子现在同时记录源比较器给出的严格事实：只要 dividend 指针仍在
  界内，当前 dividend 项次数严格小于所选 root 次数。对 dividend-source
  分支，canonical chain 的成员系数定理直接证明机器 word 转入 `ZMod`
  后等于 `toPoly.coeff`。对 heap-source 分支，证明遍历任意真实
  `dividend.toList` 成员并恢复其数组索引：指针前的成员由
  `ConsumedDividendAbove + decrease` 严格高于 root；指针处及之后的成员
  由 selector 比较与 canonical index order 严格低于 root。因此该次数
  完全不存在于 dividend，`listSum` 系数确为零。这排除了把 heap 分支的
  零初值当作无条件规格默认值。下一步把该等式与 equal-degree trace
  求和合并成 outer iteration 当前次数的 residual coefficient 等式。
- 单次 outer iteration 的 residual coefficient 等式已闭合。新增组件分解
  定理从一次成功的真实 iteration 中同时恢复同一次
  `ConsumeEqualDegree`、`Emit`、reverse `ReinsertLin` 执行及最终状态等式，
  防止只证明一个与主执行无关的辅助 consume。随后对该真实 consume
  应用 equal-degree 求和语义，并用 selector 的 dividend coefficient 桥
  改写初值，得到
  `consumed.coefficient = dividend.coeff(frontier) - Σ(qTail·dTail)`。
  乘积列表仍逐项携带真实 quotient/divisor 存储成员和次数和等于
  frontier 的证据。下一步必须加强 heap/owner 表示为卷积覆盖不变量，
  证明这份执行 trace 不仅 sound，而且对当前 quotient×divisor 尾项
  complete 且无重复。
- 卷积覆盖证明的节点位置基础已闭合。对一次真实 `ConsumeNode`，现精确
  证明：advance 分支把当前节点加入 `lin`，exhausted 分支利用源断言
  `nodeIndex = reset_h` 把它纳入扩张后的 reset 前缀；同时旧 `lin` 集合
  与旧 `reset_h` 上界都单调保持。该性质沿实际 `next` 链按 owner card
  良基递归提升：exact owner 中每个节点在完整消费结束后必定位于最终
  `lin` 或最终 reset 前缀。递归中使用真实 node update 后的 tail chain
  ownership 与实际返回 `next`，没有用节点计数相等替代成员覆盖。
  下一步把此“已消费 owner 重分类”与未消费 heap bucket 的 disjoint
  ownership 合并，建立全 node block 的 `heap ∪ lin ∪ reset` 完整分区。
- 全 node block 覆盖谓词已建立并完成入口桥。
  `PairVecDivVHCNodesCovered` 明确要求每个已分配节点落在 reset 前缀、
  `lin` 栈或某个真实 heap slot 的 exact owner 中；区域间唯一性继续由
  已有 `ResetReady + LinReady + HeapChainsOwnedAway` 提供。初始化时
  `nodes.size = divisor.size - 1 = reset_h`，因此全部节点真实位于 reset
  前缀。反向查询定理进一步证明：在覆盖、reset-ready 条件下，任何
  `mono ≠ none` 且不在 `lin` 的节点必由实际 heap bucket 持有，reset
  情形因其真实 `mono = none` 直接矛盾。此外已证明 ConsumeNode 及完整
  ConsumeChain 始终保持 node-array size，后续覆盖传递不会偷换节点域。
  下一步补强 `VHC_extract`：不仅新 head 来自旧 heap，还要证明每个非 root
  旧 head 都存活；然后即可组合 root owner 重分类得到覆盖保持。
- `VHC_extract` 的双向 head 守恒及 equal-degree 全覆盖保持已闭合。
  heap slot 唯一性先被转成 `heap.toList.Nodup`；结合 extract 精确 size
  减一、新 head 单向来自旧 heap、以及新 heap 排除唯一 root，证明新
  head finset 精确等于旧 head finset 删除 root。因此每个非 root 旧 head
  都恢复出新 heap 中的具体 slot。利用该反向存活定理，root owner 的
  所有节点由 chain 重分类定理进入最终 `lin/reset`，旧 `lin/reset` 由
  单调性保持，其他 owner 则随非 root head 留在新 heap；node-array size
  保持用于统一量化域。最后按 checked extract 的实际 heap size 强归纳，
  `ConsumeEqualDegree` 同时保持 exact ownership 与
  `PairVecDivVHCNodesCovered`。下一步把覆盖继续穿过 emit activation 和
  reverse reinsertion，使下一轮 outer loop 重新得到 `lin = []` 的完整
  heap/reset 分区。
- 覆盖不变量已提升为单一见证的组合状态。新增
  `PairVecDivVHCStateCovered`，要求同一个 `owners` 映射同时证明真实 heap
  chain ownership 和整个 node block 的 `heap/lin/reset` 覆盖；这排除了在
  insertion 改换 bucket head 后，分别挑选两个互不相关的存在见证并误称
  覆盖已传递的漏洞。初始化由空 heap 与真实 `pairVecDivVHCInit` 大小直接
  建立该状态，equal-degree 消费则复用刚闭合的 ownership+coverage 联合
  定理保留同一个见证。下一步必须让 activation/reinsert 的每次真实
  `VHC_insert` 构造新的共同见证，再沿各自的良基循环归纳。
- exact chain 的 owner 唯一性已闭合。对固定 node array 与 chain head，
  `PairVecDivVHCChainOwns` 的 owner finset 现在按真实 `next` 链良基递归证明
  唯一；由此得到任意两个合法 heap ownership 映射在每个实际 heap head
  上必相等。`NodesCovered.congr_owners` 进一步证明覆盖可以安全地换到同一
  heap/node 状态上的任意 ownership 见证，`StateCovered.covered_with` 将该
  事实封装为后续 insertion 证明的消去规则。这一步没有假定 owner 映射
  全局函数相等，只在真实 heap slot 对应的 head 上使用 chain 唯一性。
  下一步追踪 insertion 对 heap head 集合的改变：普通 bubble 保留旧 heads
  并加入新 head，merge 则以新 head 的 `next` chain 取代被合并旧 head。
- insertion 的 bubble 路径已补齐反向 head 守恒。新增通用定理证明：若
  target heap 与 source heap 等长、每个 target 值来自 source，且两边 slot
  值均唯一，则两者的 head finset 精确相等，所以每个 source head 都能恢复
  出 target 中的具体 slot。该结论分别应用于真实 `VHC_bubble` 与
  `VHC_bubble_below` 执行；证明使用两段生成循环已有的 size、values-from
  与 uniqueness 定理，不把单向复制关系误当作排列。下一步用这些反向
  存活见证传递普通插入分支的 node coverage，并单独处理 merge 分支中
  被旧 head 替换为 `newNode -> oldHead` 的 chain 扩张。
- reverse `lin` reinsertion 的完整组合覆盖已闭合。首先证明 `lin.pop`
  精确保留除末元素外的所有成员，并证明真实 `VHC_set_next` 以及所有
  `VHC_insert` 成功分支不改变 node-array size。随后将 insertion 分成两类
  精确处理：equal-degree merge 把具体 slot 的旧 head 替成 `newNode`，其
  owner 严格为 `insert newNode oldOwner`；empty/root-bubble/bubble-below
  路径则建立 singleton 新 owner，并用反向 head 守恒搬运每个旧 bucket。
  组合定理覆盖了 `VHC_insert` 的全部真实控制流。最后沿实际 `lin.pop`
  的 size 递减做良基递归，证明完整 `pairVecDivVHCReinsertLin` 返回时
  `lin = #[]` 且同一 witness 同时满足 exact heap ownership 与全 node
  coverage。下一步对 reset activation 建立对应的“reset 前缀减一并把该
  node 激活插入 heap”组合覆盖循环，再接入 outer iteration。
- reset activation 的组合覆盖循环已闭合。单步证明把真实
  `nodeIndex = resetH - 1` 从 inactive reset 前缀重分类到
  `lin.push nodeIndex`，执行生成的 `pairVecDivVHCActivate` 后证明旧 heap
  chains 不受该 fresh node 写入影响，并建立临时 lin 的 active/nodup 与
  heap-away 条件；随后直接调用已覆盖全部控制流的真实 `VHC_insert`
  定理，得到 node 进入 heap 且 reset 前缀精确缩为 `resetH - 1`。完整
  `pairVecDivVHCActivateReset` 再按实际 `resetH` 良基递归，同时传递
  away、lin-ready、reset-ready 和同一 witness 的 ownership+coverage，
  返回时 reset 前缀严格为零。这里临时 `lin.push` 只是节点位置不变量的
  表达，执行仍是原始 activate 后直接 insert，没有新增算法步骤或 fuel。
  下一步把 equal-degree consume、emit activation、reverse reinsertion 的
  三段组合覆盖接入一次完整 outer iteration，并恢复下一轮入口的
  `lin = #[] / resetH` 全分区。
- 一次完整 outer iteration 的组合覆盖已接通。`Emit` 的新定理逐分支
  对齐真实代码：仅在 nonzero coefficient、degree guard、nonzero quotient
  cell 三个条件同时成立时调用实际 `ActivateReset`，其余路径保持原
  state/resetH；activation 的 reset-ready 参数被正确保留为 append 前的
  `quotient.size`，没有错误等同于新 quotient 的 size。outer theorem 从
  `pairVecDivVHCOuterIteration_components` 恢复同一次 consume、emit、
  reinsert 执行，依次应用 equal-degree coverage、emit coverage 和完整
  reverse-lin coverage，最终证明返回状态重新满足 `lin = #[]` 的同一
  witness ownership+全 node coverage。所需 reset-ready 来自同一次 consume
  的 node invariant 定理，未用独立假设替换中间执行状态。下一步把该
  不变量沿 outer loop 的严格 frontier degree 递减传播，并利用全覆盖+
  disjoint ownership 证明每个 quotient×divisor cursor 对被处理恰好一次，
  闭合卷积 trace 的 complete/no-duplicate 方向。
- outer-loop 归纳现在强制携带组合覆盖。现有 canonical-success 定理新增
  `PairVecDivVHCStateCovered heap nodes #[] resetH` 前提，并在每次真实
  iteration 后调用完整覆盖定理，把返回的共同 owner/coverage witness
  传入严格 frontier-degree 递减的递归调用。因此后续 general-division
  refinement 无法只使用 canonical/denotes 等单向不变量而绕开 node block
  totality。重新审计顶层接口后确认：当前文件仍没有 general
  `pairVecDivIR_refines`，只有 single-divisor refinement、general executor、
  canonical preservation 与单次 residual soundness；所以 heap 基础设施的
  高完成度不能视作整个 SQF refinement 已接近闭合。下一步必须新增 cursor
  processed-prefix 不变量，把 node coverage 提升为每个 quotient×divisor
  product 在对应 frontier 恰好被消费一次，再建立 general division 的
  多项式等式，之后才能连接 SQF。
- cursor processed-prefix 不变量已建立并闭合单节点消费。
  `PairVecDivVHCCursorPrefixAbove degreeLimit` 对每个实际 node cursor 明确
  量化所有 `q < node.quotientIndex` 的真实 quotient/divisor 存储项，要求
  其乘积次数不低于当前 outer frontier bound；初始化时所有 cursor 的
  quotientIndex 均为零，因此该前缀严格为空。通用 `set_advance` 定理证明
  将 cursor 加一只新增旧 current cell，其余 prefix 与其他 nodes 通过
  concrete array set 保持。对真实 `ConsumeNode` 的 advance/exhausted 两个
  成功分支，利用 `NodeDenotes` 与当前 mono 次数精确证明新增 cell 的乘积
  次数就是正在消费的 frontier，故 processed-prefix 保持；错误分支没有
  伪造返回值。下一步沿 exact `next` chain 和 equal-degree heap loop 提升
  该不变量，再结合 canonical quotient 的 cursor 后缀严格递减证明当前
  frontier 的所有乘积只能落在当前 cursor heads 中。
- processed-prefix 已沿完整 bucket chain 提升。新定理同时携带
  `ChainAtDegree`、全局 `NodeDenotes` 与 prefix invariant，按实际
  `unvisited.erase nodeIndex` 的 finset card 严格递减：每步从 concrete
  chain 节点恢复当前 mono 次数，调用单节点 advance/exhausted 保持，使用
  `ConsumeNode_get_ne` 把 tail 的 degree 证据搬到更新数组，并用真实
  `ConsumeNode_preserves_denotes` 为递归状态重建语义。root-bucket 包装再把
  source 使用的 `Finset.range nodes.size` chain 直接接入这一归纳。下一步
  结合 exact owner subset-range 与 heap-chain homogeneity，把该结论沿每次
  checked extract 的 heap-size 严格下降提升到完整 equal-degree loop。
- processed-prefix 已沿完整 equal-degree heap loop 提升。新定理不是仅
  对 bucket 做列表归纳，而是展开真实 `ConsumeEqualDegree` 控制流：从
  root 的 exact owner 和 heap-chain homogeneity 恢复整条 source chain
  的目标次数，执行真实 root-bucket consume 后传递 prefix 与
  `NodeDenotes`，再经过 checked extract 传递 ownership/homogeneity，
  最后以返回 heap 的严格 size 下降作良基归纳。次数不匹配和
  empty heap 分支均从实际返回值直接保持不变式，没有 fuel、
  fallback 或语义神谕。下一步处理 emit 导致的 quotient append 与
  frontier bound 变化，再把 prefix 传过 activation/reinsertion 和完整 outer
  iteration。
- emit 前的 quotient append 边界已语义化。新增
  `PairVecDivVHCCursorIndicesBounded` 表达每个真实 cursor 最多处在
  quotient 的 one-past-end；从同一 witness 的 state coverage 出发，
  reset prefix 通过 `ResetReady` 得到精确等于 size，lin 与 heap-owned
  nodes 则分别通过 `LinReady`/exact chain ownership 恢复 active mono，
  再由 `NodeDenotes` 的真实 quotient lookup 得到严格小于 size。
  在该边界上已证明 quotient `push` 保持旧 processed-prefix：任何
  `q < cursor` 都不可能是新尾下标，所以 lookup 必然回到原
  quotient cell。同时加入 frontier bound 单调性，为下一轮严格
  更小 frontier 复用前缀证据。下一步证明 activate/set-next/heap
  insertion 不改变 cursor fields，把 append 后的 prefix 传过完整
  `ActivateReset` 与 reverse reinsertion。
- processed-prefix 已穿过 emit 的全部真实状态变换。通用
  `set_fields` 引理仅允许更新节点的 `mono/next`，并强制
  quotient/divisor cursor 字段与原节点相等；它已分别绑定到
  真实 `VHC_activate` 和 `VHC_set_next`。`VHC_insert` 的结论由已有
  全分支 `Insert_nodes_result` 恢复本次实际 set-next，不重建抽象
  heap 算法。完整 reverse-`lin` 循环按真实 `lin.pop` 大小良基
  递归传递 prefix；完整 `ActivateReset` 则按真实 `resetH - 1`
  依次执行 activate、insert 并递归。最后 `Emit` 的 coefficient=0、
  degree guard、value=0 均从实际 unchanged 返回保持，唯一 append
  分支先用 cursor bound 排除新尾进入旧 prefix，再调用完整
  reset activation 定理。下一步从 outer-iteration components 恢复同一
  consume/emit/reinsert 执行，组合 state-derived cursor bound 并把 prefix
  传到严格更小的下一 frontier。
- processed-prefix 已接入完整 outer iteration。新定理从同一次
  `OuterIteration_components` 恢复实际 consume、emit 和 reverse-reinsert
  结果；先用选中 frontier 对上轮的 strict upper bound 做单调
  收紧，再依赖 exact heap ownership+homogeneity 把 prefix 穿过完整
  equal-degree consume。consume 后的 cursor bound 不是作为新假设，而是
  从同一返回状态的 `StateCovered`、`LinReady`、`ResetReady` 和
  `NodeDenotes` 现场推导，然后才用于 quotient append/emit，最后
  经真实 reverse reinsertion 返回 `result.nodes/result.quotient`。返回
  prefix 的 bound 正是本轮实际 `frontier.degree`，因此可直接作为
  良基 outer-loop 下一轮的 strict upper bound。下一步将该不变式
  线程到完整 outer-loop，然后用 prefix+严格下降的 quotient 后缀
  闭合当前 frontier 卷积项的 completeness/no-duplicate。
- 当前 frontier 卷积项已完成 cursor 定位的核心排他证明。对
  任意真实 node 与其固定 divisor tail cell：若 quotient 下标严格
  早于当前 cursor，processed-prefix 与上一轮 strict frontier bound 给出
  乘积次数严格高于当前 frontier；若下标严格晚于 cursor，
  quotient canonical 的严格次数递减与 `NodeDenotes` 给出乘积次数
  严格低于当前 node mono/frontier。两侧合并后，新定理
  `pairVecDivVHCProductAtFrontier_eq_cursor` 证明任何次数恰等于
  frontier 的 quotient×divisor-tail 存储乘积，其 quotient 下标必然
  精确等于该 node 的当前 cursor。这已排除同一 divisor row 中
  prefix/suffix 遗漏或重复的可能。同时确认 outer-loop 还需显式
  重建 reinsertion 后的 heap-chain homogeneity，不会将它作为未证假设
  塞入递归。下一步用 node total coverage+fixed divisor indices 把该
  row-local 定位提升到所有 divisor-tail rows，并补齐重插入后的
  homogeneity 恢复。
- divisor-tail row 与 source node 的总映射/单射已形式化。在
  `nodes.size = divisor.size - 1` 与永久 `divisorIndex = i + 1`
  不变式下，每个 `0 < d < divisor.size` 都由真实数组槽
  `d - 1` 表示，且任何两个声称表示同一 `d` 的节点槽必然
  相等。对 total coverage 的 reset 分支，`ResetReady` 给出该 row
  cursor 精确等于 `quotient.size`；因而任意实际 quotient lookup
  的下标都严格位于 cursor 之前，由 processed-prefix+上轮 strict
  bound 得到其与该 divisor cell 的乘积次数严格高于当前
  frontier。所以 inactive/reset rows 不可能包含本轮应消费的
  frontier product，该结论来自真实 one-past-end cursor，不是忽略
  inactive nodes。下一步处理 coverage 的 heap-owned 分支：需用
  heap root dominance+chain homogeneity 证明当前 cursor mono 不高于选中
  frontier，再同 row-local 定位合并。
- heap-owned 节点到选中 frontier 的次数上界已从真实堆序推出。新增
  `pairVecDivVHCParent_lt` 证明每个非根槽的父槽严格更小，并以槽下标的
  强归纳将 `PairVecDivVHCHeapOrdered` 沿实际 parent 链提升为“任意有效槽
  次数不大于 root 次数”。随后直接展开真实 `pairVecDivVHCSelectFrontier`
  的 dividend/heap 分支，证明 root 次数不大于实际返回 frontier；组合后
  得到任意有效 heap slot 的 active mono 次数不大于本轮 frontier。
  这些结论要求 concrete heap pointer validity、heap order、成功 mono read
  和成功 selector run，没有引入 L2 oracle、fuel 或默认结果。下一步把
  state coverage 恢复出的 owner chain、chain homogeneity 和这个 slot 上界
  合并，排除 heap-owned row 的 cursor 前后两侧，从而证明每个 frontier
  product 都确实位于本轮被消费的唯一 cursor。
- heap-owned row 的 cursor 定位已闭合。新良基引理同步递归 exact
  `ChainOwns` 与 `ChainAtDegree`，证明 owner 集中的每个真实节点均 active
  且次数等于 bucket head；再用上一阶段的 slot-to-frontier 上界得到该节点
  mono 次数不大于实际 selector frontier。推广后的 row-local 排他证明只
  需要这个上界：cursor 之前的 product 由 processed-prefix 严格高于
  frontier，cursor 之后的 product 由 quotient canonicality 严格低于
  node mono、进而不高于 frontier，所以次数等于 frontier 的存储乘积只能
  精确落在当前 cursor。下一步用 `StateCovered` 在空 `lin` 状态下对每个
  divisor-tail row 分解 reset/heap-owned 两种情况，得到全 row completeness，
  再对接本轮 equal-degree consume 的 exact owner chain。
- 全 divisor-tail row 的 frontier 定位已由 `StateCovered` 闭合。在真实 outer
  边界 `lin = #[]` 下，固定的 row 槽 `d - 1` 不可能处于临时栈；若处于
  reset prefix，则此前的 reset-row 定理推出该具体 quotient/divisor cell
  乘积次数严格高于 frontier，与目标次数相等矛盾；因此它必由某个实际
  heap slot 的 exact owner chain 持有。结合 heap-owned row 定位，新定理
  返回 concrete slot/head/node、owner membership、永久 divisor index 和
  `q = node.quotientIndex`。这覆盖所有 `0 < d < divisor.size` 的实际存储
  cell，不跳过 inactive rows。下一步增强 execution-indexed `ConsumeTrace`
  的 products completeness（现有定理只证明每个 trace product 的
  soundness），再沿真实 heap-size 良基 equal-degree loop 汇总所有这些
  唯一 cursor。
- concrete bucket trace 的 products completeness 已建立。新定理直接对
  execution-indexed `PairVecDivVHCConsumeTrace` 归纳，并与 exact
  `ChainOwns` 同步：head step 的实际 quotient/divisor 数组读取就是 trace
  新增的 coefficient pair；tail step 使用真实 `ConsumeNode_next` 和
  非当前槽保持定理，把剩余 owner chain 与 `NodeDenotes` 搬到更新后的
  nodes 后递归。因此 owner 中每个节点的当前存储乘积都确实出现在 trace
  products 中，而不只是证明 trace 已有元素合法。root-bucket 包装现在
  同时返回 coefficient 等式、products soundness 和 owner completeness。
  下一步沿 checked extract 后严格减小的 heap size 汇总各 root trace，证明
  full equal-degree products 对初始所有 frontier buckets 完备。
- heap pointer validity 已不再作为 frontier completeness 的独立外部假设。
  新定理从 concrete `HeapChainOwnership` 本身推出每个实际 heap slot 都指向
  owner chain 的 active head：exact `ChainOwns` 给出 head membership，随后
  owner-member active 定理恢复真实 node/mono lookup。所有 heap-owned row
  与全 row frontier 定位定理现只要求 ownership、homogeneity 和 heap order，
  避免后续递归重复携带一个可推导的不变量。下一步仍需证明真正关键的
  `VHC_extract` heap-order preservation；这不能由 pointer validity 代替，
  也不会假设 checked extract 自动维持堆序。
- extract 堆序证明的祖先侧已闭合。新增 one-edge 定理把 concrete
  `HeapOrdered` 的 array/map 形式还原为 child mono 次数不大于 parent mono；
  随后按实际 `pairVecDivVHCParent` 严格下降做强归纳，证明任意非根 slot
  的次数都不大于 root 的某个直接子节点（槽 1 或 2）。父槽为 0 时通过
  `(slot - 1) / 2 = 0` 的自然数除法边界精确推出 slot 只能是 1/2，未把
  binary-heap 祖先关系当作算术自动事实。下一步证明 sift-down 递归只写
  当前下降路径、不会回写更早槽，并据此识别 extract 新 root 为两直接
  子节点/last sentinel 的最大者。
- sift-down 的 root 写入语义已闭合。良基 `get_before` 定理逐分支展开真实
  `pairVecDivVHCSiftDown`，证明递归只会写当前节点及其严格后代，所有早于
  当前节点的 array lookup 保持；度量仍是源码的 `limit - child`。在初始
  `i=0, child=1` 且左右子存在时，随后按源码比较分别处理 left/right
  selection 与 selected/last comparison：若递归下降，`get_before` 保证
  slot 0 保持第一次写入；若停止，slot 0 就是 last sentinel。所得真实
  root head/mono 同时支配 left、right、last 三个候选。下一步把 sift result
  的 slot 0 穿过真实 `pop`，并结合所有旧非根 slot 到 left/right 的上界，
  得到 checked extract 后 root 对全部 surviving heads 的支配。
- `pop`/extract 层的新 root 全局支配已闭合。长度二的真实执行单独展开，
  证明 sift child guard 直接把 last 写到 root 后 pop；长度至少三时将上一
  阶段的 left/right/last 支配穿过 `Array.pop`。对任意返回 slot，真实
  `Extract_valuesFrom` 找回旧 slot，root-exclusion 与 heap-head uniqueness
  排除旧 slot 0，再用旧 slot 到直接 root child 的支配和新 root 对两个
  child 的支配传递得到最终上界。新定理因此证明成功 extract 的新 root
  支配全部 surviving heads，而不仅是输出值来自旧 heap。下一步建立实际
  consume/extract 桥：旧 heap order 在消费前 nodes 上成立，sift 在消费后
  nodes 上运行；需用 root-owner 与其他 owner 的 disjoint 更新保持证明，
  将所有非根 head mono lookup 搬到消费后 nodes 后应用本定理。
- consume→extract 的双 storage 桥已闭合。`Extract_root_dominates` 现显式
  区分提供旧 heap order 的 `sourceNodes` 与实际执行 sift 的 `nodes`，只
  要每个非根 heap head 的 mono 成功读取在两者间等价即可。该等价由真实
  root-bucket consume 推出：root owner 与其他 owner disjoint，consume
  trace 的 unvisited 集覆盖其他 owner，且每个 unvisited node lookup
  原样保持；特别得到所有非根 head 的 mono iff。组合定理现直接证明
  `ConsumeRootBucket` 后以 `bucket.nodes` 执行 `VHC_extract`，新 root 支配
  返回 heap 的全部槽，没有要求已消费旧 root 仍 active。下一步需要将
  这一 root-max 结果加强为完整 parent/child heap order preservation；否则
  第二次及后续 extract 仍无法合法复用 heap-order 前提，不能提前宣称
  full equal-degree completeness。
- sift 的候选支配已从初始 root 推广到每个递归 hole。新定理接受实际
  `i/child/limit/lastNode`，在 `child < limit < heap.size` 下逐分支展开同一
  source recursion；selected child 继续下降时，早先槽保持定理确保更深
  recursion 不会回写当前 `i`，停止时则直接读取 `heap.set i lastNode`。
  因而每层最终写入 hole 的 mono 同时支配该层 left、right 与 saved last。
  这是完整 heap-order preservation 的局部核心：下一步对递归路径上的
  parent/child 边使用此定理，对路径外边使用 concrete set/get 保持，按
  `limit - child` 归纳组合成 sift 后全部父子边有序。
- heap order 已建立可逆的 degree/read 前缀表示。
  `PairVecDivVHCHeapDegreesOrderedUpTo limit` 对每个 `child < limit` 直接量化
  concrete child/parent heap head lookup 与两次成功 mono read，表达同一
  次数不等式；旧 `HeapOrdered` 在 `limit ≤ heap.size` 时推出它，反向在
  `limit = heap.size` 时通过逐个还原 node lookup/active mono 与 map/join
  证据无损重建原谓词。该转换不是降低验证标准，而是把 sift 中的
  `Array.set` 局部效应从嵌套 Option map 形式解耦出来。下一步在这个等价
  前缀谓词上证明单次 child-copy 与最终 last-write 保持所有未受影响边，
  再转回现有 `HeapOrdered` 供 equal-degree 递归使用。
- 单槽更新的 heap-order frame 与 leaf write 已闭合。通用
  `set_of_affected` 把 `(heap.set i newHead)` 的证明精确分成
  `child = i`、`parent child = i` 两类受影响边；其他边通过两个真实
  `Array.getElem?_set_ne` lookup 还原到旧前缀堆序。`set_leaf` 再处理没有
  active children 的最终 write：唯一可能受影响的是 `i → parent i`，从
  新 mono 到旧 parent 的显式上界恢复；若有边声称 parent 为 leaf，则与
  no-children 证据矛盾。下一步实例化非 leaf 的 child-copy：用 selected
  对左右候选的支配证明两条下边，用旧 child≤parent 链证明写入值仍不
  超过 hole 的旧 parent，然后递归到 selected hole。
- 非 leaf 单槽写入的通用父/子 frame 已闭合。首先从源码 parent 公式
  `(child - 1) / 2` 精确证明其逆像只有 `2*i+1` 与 `2*i+2`，包含自然数
  除法上下界证明。随后 `set_parent` 将受影响边分开：`child=i` 的向上边
  由新 mono 的显式 parent 上界恢复；`parent child=i` 的向下边先把
  parent lookup 对齐到新 head/mono，再将 child lookup 通过 set-ne 还原，
  使用新 mono 对该 child 的显式支配。结合上一阶段 frame，所有其他边
  自动保持。下一步用 source selection 的 left/right 比较实例化 `hdown`，
  并从旧堆序的 `selected→i→parent(i)` 链实例化 `hup`，得到实际
  child-copy 后的前缀堆序。
- 真实 selected-child copy 的前缀堆序保持已闭合。引理显式接收 hole、
  left、right 和 selected 的 heap head 及各自真实 lookup/
  `pairVecDivVHCMono` 结果；selected 必须是左右子之一，且比较分支
  必须给出它对两个候选的 degree 支配。向上边由旧堆序的
  `selected → hole → parent(hole)` 传递恢复，向下边用 parent 逆像
  只能是 left/right 的定理逐支恢复，未受影响边继续由 set frame
  保持。这一步使用的是实际 `Array.set i selectedHead`，没有以 L2
  结果回填或使用 fuel。下一步把该单步与 leaf last-write 沿
  `pairVecDivVHCSiftDown` 的良基递归组合，得到整段 sift 的堆序保持。
- 整段生成的 `pairVecDivVHCSiftDown` 已证明保持 active-prefix degree
  heap order。证明与实现使用同一个 `termination_by limit - child`
  的良基度量，逐个展开 left/right selection 和 greater/stop 分支。
  递归分支先用真实 selected-child `Array.set` 保持堆序，再把
  `lastMono < selectedMono` 作为下一个 hole 的向上不变式；stop
  分支则直接证明 `Array.set i lastNode` 的上下边。特别单独
  处理了 `right = limit`：生成代码确实读取该 sentinel，但它不被
  偷算为 active heap child，只有真正选中并递归的节点才要求
  `< limit`。下一步从 extract 入口的旧堆序和 last-slot lookup 构造
  该定理的初始 root/parent 前提，然后经 `pop` 转回完整
  `PairVecDivVHCHeapOrdered`。
- `VHC_extract` 的完整 heap-order preservation 已闭合。从 heap pointer
  validity 取出真实 root/last node 及 active mono，用旧堆序的
  parent 链证明 saved last 不超过旧 root，以此初始化上一阶段
  良基 sift 定理。得到 shifted 的 active-prefix 堆序后，新增
  `HeapDegreesOrderedUpTo.pop` 通过真实 `Array.getElem?_pop` 将每条
  child/parent lookup 还原到 shifted，最后无损转回现有
  `PairVecDivVHCHeapOrdered`。单元素 heap 的空结果也由真实 extract
  size 定理单独覆盖。因此后续每次 extract 都可重新获得完整
  堆序前提；下一步将它与 consume 后 nodes 的非根 mono 保持桥
  组合，恢复重复 equal-degree consume 的不变式。
- consume 后的非根 heap edge 堆序已从消费前 nodes 严格搬运到
  `bucket.nodes`。证明不假设已消费 root 仍 pointer-valid；它只对
  `parent child > 0` 的真实 heap lookup 使用 root-owner 与其他 owner
  disjoint 推出的 mono iff，将 child/parent 两个成功读取还原到
  消费前 nodes，再应用旧完整堆序。这为 consume→extract 组合
  建立了正确分解：新 root 边由 extract root-dominance 处理，非根边
  由本定理和 sift path frame 处理，而不会把失活 root 塞回不变式。
- consume→extract 与完整 equal-degree 良基循环的 heap order 已闭合。
  新的 root-sift 定理不要求已消费旧 root 仍 active：第一次真实
  `Array.set 0 selected/last` 直接由 left/right/last 比较重建两条
  root 边，非根边由上阶段 consume 搬运定理保持；如果继续
  下沉，新 hole 已是 active selected child，因而无损进入完整良基
  sift 定理。该结果经真实 pop 转回 `HeapOrdered`，再与 exact
  chain ownership 一起沿 `heap.size` 严格下降的
  `pairVecDivVHCConsumeEqualDegree` 归纳，使每一次重复 extract 都拥有
  真实完整堆序前提。下一步用这个循环不变式将已有的
  root-dominance 与 bucket product completeness 提升为整段 equal-degree products
  completeness。
- 整段 equal-degree products completeness 已闭合。新谓词
  `PairVecDivVHCProductsCoverDegreeOwners` 对每个初始 heap slot/head、
  每个处于目标 degree 的 owner chain 及其每个具体 node，要求该
  node 真实 cursor 指向的 quotient/divisor 系数对出现在 products 中。
  强定理沿实际 `heap.size` 良基下降：当 root degree 匹配时，
  root owner 由 bucket trace 的逐节点 coverage 处理，其他 owner 通过
  disjoint-chain 的精确 lookup 保持和 extract surviving-head 见证进入递归；
  当 root degree 不匹配时，完整 heap order 与全局 degree 上界排除
  任何更深的目标-degree chain。最终同一个 `rootProducts ++
  tailProducts` 同时满足模系数等式、每项 soundness 和所有 owner
  行的 completeness，不是两个无关的存在性列表。下一步将
  frontier-row owned-cursor 定理与该 coverage 组合，将 products value 替换为
  L2 polynomial multiplication 在 frontier degree 的系数。
- products coverage 已从集合式 membership 进一步加强到精确重数语义。
  新定义 `pairVecDivVHCNodeProductValue` 直接从真实 node cursor 读取
  quotient/divisor 系数并在 `ZMod p` 中相乘。对任意具体
  `PairVecDivVHCConsumeTrace`，新定理证明 `ProductsValue products`
  精确等于 owner Finset 上的 node-product 求和。证明沿 trace 归纳，
  用 `owner = insert nodeIndex (owner.erase nodeIndex)` 分离当前行，并用
  consume-node 对其余 owner lookup 的精确保持搬运尾和。因此即使
  多个不同行恰好产生相同系数对，也会按行数重复计入，不会把
  membership coverage 误当成 multiset equality。下一步沿 equal-degree
  extract 循环对这个 owner-sum 等式做 disjoint union，再与所有 frontier
  rows 的唯一 cursor 对应连接。
- extract 对 heap-owned 节点有限集的精确差分已闭合。新定义
  `PairVecDivVHCHeapOwnedNodes` 将真实 `heap.toList` 中每个 head 的
  owner Finset 取 `biUnion`。对任意成功 `VHC_extract`，证明返回
  heap 的 owner 并集精确等于旧并集减去 `owners heap[0]`。正向使用
  extract values-from 和 root exclusion，反向使用每个非根 head 的具体
  surviving slot 见证，所以这不是只能支持 soundness 的子集关系。
  配合 owner 之间的 disjointness，下一步可将每个 root trace 的精确
  owner-sum 沿 equal-degree 循环分解为不交求和。
- 目标-degree heap-owned 节点集的精确递推已闭合。新的
  `PairVecDivVHCNodeAtDegree` 以真实 `pairVecDivVHCMono` 成功读取
  定义节点 degree，`HeapOwnedNodesAtDegree` 则在真实 owner 并集上
  过滤。consume→extract 后的该集合精确等于旧集合减去 root owner：
  非根 owner 的节点 lookup 通过 unvisited/disjoint-chain 证明前后完全
  相同，所以 degree 过滤也前后等价。另外，当真实 root mono
  degree 命中目标时，chain homogeneity 与 exact ownership 证明整个
  `owners heap[0]` 都包含于目标集。下一步可直接把 root trace 的
  owner-sum 与递归尾的 set-difference sum 合并为初始目标集的完整求和。
- equal-degree 循环的非匹配停止分支已证明目标 owner 集为空。
  对任意假设存在的目标-degree owned node，exact ownership 找到其
  heap head，chain homogeneity 将该 node degree 与 head degree 对齐，完整
  heap order 再给出 `head ≤ root`。结合所有 head 的全局 `≤ target`
  上界，得到 root degree 必等于 target，与实际循环的
  `rootMono.deg ≠ degree` 停止 guard 矛盾。这闭合了全循环 owner-sum
  归纳的 base/stop 情形，不需要任何规格层“无剩余乘积”假设。
- 整段 equal-degree 循环的精确重数 owner-sum 已闭合。
  `pairVecDivVHCConsumeEqualDegree_products_complete` 返回的同一个
  products 列表现在同时证明：执行系数等式、每个产品的
  degree soundness、每个 owner 行的 coverage，以及 `ProductsValue`
  精确等于所有初始目标-degree owned nodes 的 node-product 有限和。
  递归分支将 root trace 的 owner sum 与尾循环的目标集差分和组合，
  用 root-owner subset 和 `Finset.sum_sdiff` 恢复初始完整和；同时用
  consume 对非根 owner lookup 的精确保持，把尾和从 `bucket.nodes`
  无损搬回初始 nodes。空 heap 与 root guard 不匹配两个 base 分支
  均已证明目标集为空。下一步只需将这个以 cursor nodes 索引的
  有限和通过 frontier-row 唯一性改写为 L2 polynomial multiplication 系数。
- 纯 L2 的稀疏列表乘积系数桥已闭合。新定义
  `pairVecDivVHCListProductCoeffValue p degree xs ys` 对 `xs × ys`
  的每个具体索引项对逐一遍历，仅在 monomial degree 之和命中目标
  时累加系数乘积。通过先对固定 quotient monomial 归纳 divisor
  rows，再归纳 quotient rows，证明该有限和精确等于
  `(listSum p xs * listSum p ys).coeff degree`。已给出完整
  quotient/divisor `toPoly` 版本和 divisor-tail 版本。该引理不依赖
  系数对值的唯一性，因而天然保留所有索引重数。下一步是唯一
  剩余的双射：目标-degree owner cursor nodes 与 `quotient.toList ×
  divisor.toList.tail` 中 degree 命中项对一一对应。
- 双射的 L2 项对→C++ owner 方向已下沉到精确有限集成员
  定理 `pairVecDivVHCFrontierProduct_mem_ownedNodesAtDegree`。对任意
  真实 quotient 索引和非首 divisor 索引，若次数和等于当前
  frontier，则固定分配行 `d - 1` 必属于生成 C++ 堆循环实际
  消费的 `HeapOwnedNodesAtDegree`；证明同时使用 nodes 全覆盖、
  reset 排除、堆所有权、行次数齐性、固定 divisor 指针和
  cursor-prefix 不变式，并保留具体节点索引，不会因相同系数值
  合并重数。下一步是证明 owner 方向的唯一项对并改写有限和。
- owner cursor nodes 与 L2 索引项对的真正双射已闭合。新的
  `PairVecDivVHCTargetPairsAtDegree` 以 `(quotientIndex,
  divisorIndex)` 而非系数值保存目标次数项对，且只包含 divisor
  tail。`pairVecDivVHCOwnedNode_sourcePairAtDegree` 从每个 owner
  节点恢复真实源索引、次数命中与相同系数积；反向定理从
  每个索引项对恢复固定行节点。
  `pairVecDivVHCHeapOwnerSum_eq_targetPairSum` 使用显式左右逆
  映射的 `Finset.sum_bij'` 证明两个有限和完全相等，因而
  同时证明无遗漏、无重复和索引重数保持。下一步只剩将该
  indexed pair sum 展开成已证明的 L2 list product coefficient。
- indexed pair sum 到 L2 polynomial coefficient 的改写已闭合。先证明
  list 双重 `foldr` 等于以 `Fin` 索引的双重和，再对 divisor
  tail 使用显式 `d ↔ d + 1` 双射，将其与 C++ 中非首
  divisor 行的 `Ico 1 divisor.size` 严格对齐。由此
  `pairVecDivVHCHeapOwnerSum_eq_productCoeffTail` 直接证明堆所有者
  消费和等于 `(quotient * divisor.tail).coeff frontier.degree`。
  新的 `pairVecDivVHCOuterIteration_residual_coefficient_toPoly` 已将这个
  结果接入真实 `pairVecDivVHCOuterIteration` 执行，得到其消费
  结果系数等于 dividend 系数减去完整 L2 tail-product 系数。
  该路径未使用规格 oracle、L2 回退或 fuel。下一步是将首项发射
  部分合并进完整 divisor 乘积不变式，并闭合外层除法等式。
- 真实 emit 分支的商多项式增量语义已闭合。
  `pairVecDivVHCEmit_toPoly_of_lead_le` 对 divisor 首项可整除
  frontier 次数的情形，证明运行后商精确增加
  `monomial (frontier.degree - leadDegree) (consumed / leadCoeff)`。
  证明展开了同一个生成 `nmod_inv_ir`/`nmod_mul_ir` 计算，
  并复用 single-term 真实精化定理证明域除法语义及非零性；
  系数为零时则证明源分支不修改商，与零单项式增量一致。
  因此 emit 的首项贡献已不再只有 canonical 属性，而是已有
  可直接合并到完整 divisor 乘积系数的 L2 等式。
- 单次外层迭代的完整乘积系数等式已闭合。先通过
  `pairVecDivVHCQuotient_mul_lead_coeff_eq_zero` 使用 quotient-above
  不变式排除旧商与 divisor 首项在当前 frontier 的贡献；
  再通过 `pairVecDivVHCEmitted_monomial_mul_tail_coeff_eq_zero`
  使用 canonical divisor 的严格降次性，排除新商单项式与
  divisor tail 在该次数的贡献。对剩余首项乘积使用有限域
  `div_mul_cancel₀`，与上一阶段的 heap tail-product 消费等式合并。
  最终 `pairVecDivVHCOuterIteration_product_coefficient` 证明真实
  `pairVecDivVHCOuterIteration` 运行后
  `(result.quotient * divisor).coeff frontier.degree = dividend.coeff
  frontier.degree`。下一步将该单步等式提升为良基外层循环
  的全局 processed-degree 不变式。
- 已证明真实外层迭代不会破坏高于当前 frontier 的乘积系数。
  `pairVecDivVHCEmitted_monomial_mul_divisor_coeff_eq_zero_above`
  将新发射商单项式分别对 divisor 首项和严格降次 tail 展开，
  证明其乘积在任意更高次数上的系数为零；由此
  `pairVecDivVHCEmit_product_coeff_above` 证明真实 emit 保持这些
  系数，`pairVecDivVHCOuterIteration_product_coeff_above` 再经由
  select/consume/reinsert 的真实执行分解，将保持性提升到完整
  单次 outer iteration。该结果与当前 frontier 的系数闭合定理
  共同构成全局 processed-degree 归纳的单步核心。下一步仍需证明
  frontier 之间可能出现的次数空隙也为零，再将不变式穿过良基
  `pairVecDivVHCOuterLoop`，最后用整除性消去低于 divisor 首项的余式。
- 已闭合稀疏 outer-loop 不跳过来源次数的核心空隙引理。
  `pairVecDivVHCTargetPairsAtDegree_eq_empty_of_gap` 不把堆替换为规格
  集合，而是对每个真实 divisor-tail 行使用 `StateCovered` 分类：
  reset 行的现有 quotient 索引严格位于 cursor 前，因此乘积次数
  至少为旧 loop bound；heap-owned 行则按索引在 cursor 前、等于
  cursor、或在 cursor 后三分，分别由 cursor-prefix、selector 最大值
  和 canonical quotient 降次性排除 `(frontier, oldLimit)` 空隙。
  由有限 indexed pair sum 桥得到
  `pairVecDivVHCTail_product_coeff_eq_zero_of_gap`。同时
  `pairVecDivVHCDividend_coeff_eq_zero_of_gap` 对 selector 的 dividend/
  heap 两种真实来源分支证明 dividend 在同一空隙的 L2 系数为零。
  下一步还需携带 quotient 首项贡献的 processed-bound 不变式，才能
  把这两条零系数结论与单步当前次数闭合、高次数保持合并为完整
  良基循环系数归纳。
- quotient 首项的 processed-bound 已作为独立机器状态不变式补齐。
  `PairVecDivVHCQuotientLeadAbove degreeLimit leadDegree quotient` 记录每个
  已发射商项与 divisor 首项相乘后的次数至少为当前旧界限。
  `pairVecDivVHCEmit_quotient_eq_or_push` 直接展开真实 emit，证明输出商
  要么不变，要么仅 push 当前 `frontier - leadDegree` 项；由此
  `pairVecDivVHCEmit_preserves_quotientLeadAbove` 与完整 outer-iteration
  版本证明良基递归下降时该不变式更新到当前 frontier。
  `pairVecDivVHCQuotient_mul_lead_coeff_eq_zero_below_processed` 排除旧商
  首项乘积落入次数空隙，最终
  `pairVecDivVHCProduct_coeff_eq_zero_of_gap` 将首项与前一提交的真实
  heap tail 结论合并，证明完整 `quotient * divisor` 在
  `(frontier, oldLimit)` 的系数为零。至此单步全次数覆盖已齐，下一步
  是用强归纳将“旧界限以上相等、空隙为零、frontier 当前相等”贯穿
  `pairVecDivVHCOuterLoop`。
- 单步全区间系数扩展已闭合。
  `PairVecDivVHCProductAgreesAbove` 表示 `quotient * divisor` 与 dividend
  在某个次数界限以上逐系数相等；
  `pairVecDivVHCOuterIteration_extends_productAgreement` 对任意目标次数
  作三分：等于 frontier 时使用真实 consume/emit 的完整系数定理，
  严格位于 frontier 与旧界限之间时使用 dividend、lead-product 与
  heap-tail 三条空隙排除定理，高于旧界限时使用旧相等不变式；真实
  emit 的高次数保持把三种情形统一到 iteration 输出商。
  因而单次生成 C++ body 已能把 L2 相等区间从 old limit 严格扩展到
  selected frontier。全局强归纳现在剩下的是证明 activation/insertion/
  reinsertion 对 heap order 与 equal-degree chain homogeneous 的保持，
  以便每个递归状态都满足该单步定理的表示前提；这些前提不会被
  当作假设跳过。
- `VHC_insert` 的节点写入与 heap 排列责任已正式解耦。
  `pairVecDivVHCSetNext_preserves_mono_read` 直接展开生成的 checked
  `set_next`，证明无论读取被写节点还是其他节点，比较器看到的
  monomial 完全不变；`pairVecDivVHCSetNext_preserves_heapOrdered` 因而
  将任意已建立的 heap order 穿过真实 `next` 写入。进一步的
  `pairVecDivVHCInsert_nodes_preserve_heapOrdered` 使用 insertion 的真实
  `nodes_result` 分解，将这个事实应用到完整 insertion 返回节点数组。
  这样剩余 heap-order 证明只需验证 `Bubble`/`BubbleBelow` 对 heap 槽位
  的父子次数关系，不再把节点链指针副作用混入数组排列证明。
- root-bubble 的终点语义已从真实递归执行中提取。
  `pairVecDivVHCBubble_stop_get` 对 generated pointer-copying recursion 作
  与源码相同的良基递归证明：成功返回时 `stop` 槽必精确包含待插入
  node。`pairVecDivVHCInsert_root_of_greater` 随后展开 insertion 的
  `newDegree > rootDegree` 分支、真实 `set_next` 与 `Bubble` 调用，证明
  返回 heap 的根确为新 node。这为 new-root 分支的 max-heap 证明固定
  了根键；下一步需证明沿祖先路径下移的旧节点保持所有局部父子边。
- upward-bubble 返回 heap 的全局根上界已闭合。
  `pairVecDivVHCBubble_new_root_bounds_all` 使用真实 bubble 的
  `ValuesFrom` 定理，将每个返回槽位反向追溯到 `oldHeap.push newNode`：
  若来自旧 heap，则由原 `HeapOrdered` 的 slot-to-root 定理及
  `oldRoot ≤ newMono` 得到上界；若来自 push 的最后槽，则通过实际
  `VHC_mono` 读取证明其键就是 `newMono`。
  `pairVecDivVHCInsert_new_root_bounds_all` 将该结果接回完整 generated
  insertion 的 greater-root 分支。现在该分支已有“返回 root 是新节点”
  及“所有返回键不超过新根”两部分；剩余局部责任是证明非根父子边
  在祖先路径复制后仍有序。
- upward-bubble 的单次祖先复制保持完整 max-heap order。
  `pairVecDivVHCSet_child_to_parent_preserves_heapOrdered` 直接针对生成实现
  的 `heap.set slot parentHead`：新槽与其父节点键相同；其子节点先由旧
  heap order 小于等于被覆盖的旧键，再由旧槽到父槽的边传递到复制的
  新键；所有未受影响边由通用 `set_parent` 框架保留。该证明使用真实
  数组读取、有效节点的 `VHC_mono` 成功语义和两级父子次数关系，没有
  引入 L2 oracle、fuel 或额外语义假设。下一步把此单步引理沿
  `pairVecDivVHCBubble` 的良基递归组合，并在根写入处使用已证明的全局
  新根上界，从而关闭 greater-root insertion 的完整 heap-order 保持。
- greater-root `Bubble` 的完整 max-heap order 保持已闭合。
  新增的 `PairVecDivVHCHeapBoundedBy` 及真实 `set`/`push` 保持引理记录
  当前 heap 中每个活动键均不超过待插入键；
  `pairVecDivVHCPush_parent_preserves_heapOrdered` 证明 append 后第一次把
  父指针复制到末槽立即恢复 heap order。
  `pairVecDivVHCBubble_to_root_preserves_heapOrdered` 随生成函数的
  `termination_by i` 做良基递归：每次祖先复制使用上一提交的局部保持
  定理，最终 `i = 0` 的真实根写入则逐条验证所有受影响根边。
  `pairVecDivVHCBubble_new_root_preserves_heapOrdered` 从原 heap order 推出
  全局新键上界并组合首次复制与递归尾部，覆盖完整 generated
  greater-root 调用。全程没有 fuel、规格 oracle 或 L2 回退。
  下一步把该定理接入 `pairVecDivVHCInsert` greater-root 分支，并处理
  `FindAnchor`/`BubbleBelow` 的非根 insertion 分支。
- greater-root heap-order 结果已接回完整 generated `VHC_insert`。
  `pairVecDivVHCInsert_newRoot_preserves_heapOrdered` 展开真实分支判定、
  `set_next newNode none` 与 `Bubble` 执行，先在写入前节点数组上应用
  完整 root-bubble 定理，再用已证的 monomial-read 不变性把 heap order
  运输到 insertion 实际返回的节点数组。由此不只 standalone bubble，
  而是完整 C++ insertion greater-root 分支已有表示不变量结论。
  下一步集中闭合 `FindAnchor`/`BubbleBelow` 非根分支。
- `FindAnchor` 到 `BubbleBelow` 的真实比较路径已闭合为 heap-order 证明。
  `PairVecDivVHCFindAnchorTrace` 精确记录 generated anchor search：每个
  climb 都带有实际 heap/node 读取及 `oldDegree < newDegree`，stop 则带
  `newDegree ≤ anchorDegree`；`pairVecDivVHCFindAnchor_trace` 从成功 raw
  执行递归地产生该 trace，而不是把路径比较当作外部前提。
  trace 的 `set_above`/`push` 引理证明后续 heap 指针写入不会伪造或丢失
  尚未消费的祖先比较。
  `pairVecDivVHCBubbleBelow_trace_preserves_heapOrdered` 沿 generated
  `BubbleBelow` 的良基递归消费 climb：祖先复制使用真实 `Array.set`
  保持定理，最终写入用 climb 的严格下界支配旧槽子树、用 stop 的上界
  保持 anchor 父边。`pairVecDivVHCBubbleBelow_push_trace_preserves_heapOrdered`
  进一步覆盖从 appended last slot 开始的完整调用，包括 anchor 就是
  first-parent 的零次复制分支。至此非根 bubble 本体不再缺 heap-order
  语义；下一步将 trace 从完整 `VHC_insert` 的 `FindAnchor` 执行接入其
  unequal-anchor 分支，并处理 equal-anchor 的链式替换。
- 完整 generated `VHC_insert` 的 max-heap order 保持已闭合。
  `pairVecDivVHCSet_sameDegree_preserves_heapOrdered` 证明 equal-degree
  bucketing 用同次数新 head 替换 root/anchor 时，上下所有父子边均保持；
  随后的 `set_next` 只改链字段，已由 monomial-read 不变性运输到实际
  返回节点数组。`pairVecDivVHCInsert_bubbleBelow_preserves_heapOrdered`
  则展开 unequal-anchor 分支，把真实 `FindAnchor` 执行生成的 trace 接到
  完整 appended `BubbleBelow` 证明。
  最终 `pairVecDivVHCInsert_preserves_heapOrdered` 对 empty、equal-root、
  greater-root、equal-anchor、unequal-anchor 逐个展开 raw 成功路径，给出
  任意成功 insertion 的统一 heap-order 定理。至此 insertion 的排列
  不变量不再是 outer-loop 的缺口；下一步证明 equal-degree `next` 链
  homogeneous 与指针有效性的完整 insertion 保持，并装配主 invariant。
- reverse-`lin` reinsertion 的完整 heap-order 保持已闭合。
  `pairVecDivVHCReinsertLin_preserves_heapOrdered` 沿 generated
  `pairVecDivVHCReinsertLin` 的 `lin.size` 良基递归：每步从实际
  `HeapChainOwnership` 推出 heap-pointer validity，调用完整
  `pairVecDivVHCInsert_preserves_heapOrdered`，再用已有 protected-set
  ownership 与 `LinReady.pop_after_insert` 为下一递归状态重建前提。
  因此 consume 后暂存到 `lin` 的节点全部按源码逆序插回时，不会在任何
  中间步骤丢失 heap order；没有执行预算或规格层替代循环。下一步把
  consume 的有序性、emit/activation 的非干扰性与此 reinsertion 定理
  合成为 outer iteration 的 heap-order 保持。
- activation/reset 的完整 heap-order 保持已闭合。
  `pairVecDivVHCActivate_preserves_mono_read_ne` 从真实 checked activation
  的数组写入证明其他节点的 comparator 读取不变；结合 ownership 对
  inactive reset 节点的 fresh-head/fresh-chain 结论，
  `pairVecDivVHCActivate_preserves_heapOrdered_of_freshHead` 保留旧 heap 的
  每条边。`pairVecDivVHCActivateInsert_preserves_heapOrdered` 随后接入完整
  insertion 定理，而 `pairVecDivVHCActivateReset_preserves_heapOrdered`
  沿 generated `resetH` 良基递归组合每一次 activate+insert。
  这关闭了 emit 的 activation 阶段可能改变比较键的真实表示风险；
  下一步直接装配完整 outer iteration 的 heap order。
- 完整 outer iteration 的 heap-order 保持已闭合。
  `pairVecDivVHCEmit_preserves_heapOrdered` 对 emit 的四条真实控制流分别
  证明：不发射分支原样保持 consumed heap；发射分支由完整
  `ActivateReset` 定理保持。`pairVecDivVHCOuterIteration_preserves_heapOrdered`
  使用 `OuterIteration_components` 提取 generated consume、emit、reinsert
  的实际成功执行，先由 equal-degree consume 得到 ownership/order 与
  reset-ready，再用 protected `lin` 不变量穿过 emit，最后调用 reverse
  reinsertion 定理得到 iteration 返回 heap 的 order。
  至此 outer loop 每轮所需的 heap 排列前提可以递归传递；下一步将它与
  homogeneous ownership、cursor prefix、processed quotient lead 和系数
  agreement 打包进全局良基递归证明。
- equal-degree bucket 的底层 homogeneous 写入语义已闭合。
  `pairVecDivVHCSetNext_chainAtDegree_insert` 直接证明真实 `set_next` 将新鲜
  节点接到既有同次数 tail 后，整个扩展链仍在指定次数；
  `pairVecDivVHCHeapChainsHomogeneous_push_fresh` 覆盖 `next := none` 的新
  单节点 bucket，并用 owner freshness 将其他 bucket 穿过节点写入；
  `pairVecDivVHCHeapChainsHomogeneous_merge_fresh` 覆盖 equal-root/anchor
  合并，显式使用新旧 head 的次数相等，而非全局假设所有节点均 active。
  `PairVecDivVHCHeapChainsHomogeneous.of_valuesFrom` 则允许 bubble 只重排
  heap head 时运输该性质。下一步将三者按完整 insertion 分支组合，再
  沿 activation/reinsertion/outer iteration 递归传递。
- 完整 generated insertion 的 ownership+homogeneous 保持已闭合。
  `pairVecDivVHCInsert_preserves_heapChainsHomogeneous_of_fresh` 与真实
  insertion 控制流逐支对应：empty/greater/unequal 分支建立 `{newNode}`
  单节点 owner，equal-root/equal-anchor 分支建立
  `insert newNode oldOwner`；bubble 分支通过真实 `Bubble`/`BubbleBelow`
  values-from 运输 head-chain 的次数性质。新节点的比较 monomial 由
  checked `VHC_mono` 与实际节点读取相等推出，equal 分支的次数相等
  直接来自源码分支条件。下一步沿 `ActivateReset` 和 reverse `lin`
  递归携带该定理，形成 outer iteration 的新 owner+homogeneous 对。
- homogeneous ownership 已穿过 activation 与 reverse reinsertion。
  `pairVecDivVHCActivate_preserves_heapChainsHomogeneous_of_fresh` 证明 inactive
  reset 节点的真实 monomial 写入对所有现有 owner chain 不可见；
  `pairVecDivVHCActivateInsert_preserves_heapChainsHomogeneous` 将其接到完整
  insertion，`pairVecDivVHCActivateReset_preserves_heapChainsHomogeneous`
  再沿 `resetH` 良基递归返回每步更新后的精确 owner map。
  对 reverse `lin`，`PairVecDivVHCHeapChainsHomogeneous.congr_owners` 利用
  `ChainOwns` 的 owner 唯一性，证明不同辅助定理产生的 owner 函数在实际
  heap head 上一致；`pairVecDivVHCReinsertLin_preserves_heapChainsHomogeneous`
  因而能沿 `lin.size` 良基递归同时传递 ownership 与 homogeneous。
  下一步将 consume、emit activation 与 reinsertion 的 homogeneous 结果
  合成 outer iteration，并进入全局 outer-loop 系数归纳。
- 完整 outer iteration 的 ownership+homogeneous 保持已闭合。
  `pairVecDivVHCEmit_preserves_heapChainsHomogeneous` 对 emit 的真实分支
  分别保留原 owner map或采用 `ActivateReset` 产生的新 map；
  `pairVecDivVHCOuterIteration_preserves_heapChainsHomogeneous` 从实际
  components 提取 consume/emit/reinsert 执行，逐阶段用 ownership 唯一性
  协调 existential owner witnesses，并最终返回 iteration 结果 heap/nodes
  上的精确 ownership 与 homogeneous。至此单轮系数桥所需的 heap order
  与 homogeneous 均可递归传递。下一步完成全局 `OuterLoop` 的良基
  `ProductAgreesAbove` 归纳及 done-state 全次数结论。
- generated outer loop 的 terminal coefficient 语义已闭合。
  `pairVecDivVHCDividend_coeff_eq_zero_of_done` 从真实 source index 已到数组
  末尾及 consumed-prefix 不变量证明当前界以下没有 dividend 项；
  `pairVecDivVHCTail_product_coeff_eq_zero_of_done` 从空 heap、`StateCovered`、
  `ResetReady` 和 cursor prefix 证明每个 divisor-tail 行都已经越过全部商项，
  因而当前界以下不存在目标乘积；
  `pairVecDivVHCProduct_coeff_eq_zero_of_done` 再与 processed quotient/lead
  不变量组合，得到完整 `quotient * divisor` 的低次系数为零。该结论不引用
  可除性规格、预期商或任何 L2 回退；下一步把它作为 `OuterLoop` 良基归纳
  的 done 分支，并在 step 分支递归传递全部表示不变量。
- VHC 节点块长度已沿完整 generated iteration 保持。
  `pairVecDivVHCInsert_nodes_size` 复用 insertion 的真实 `SetNext` 执行见证，
  证明所有 heap 分支都只改写预分配节点；随后
  `pairVecDivVHCConsumeEqualDegree_nodes_size` 沿 heap-size 良基递归组合
  bucket consume，`pairVecDivVHCActivateReset_nodes_size` 与
  `pairVecDivVHCReinsertLin_nodes_size` 分别沿 reset/lin 的真实良基递归组合
  activate+insert，`pairVecDivVHCEmit_nodes_size` 覆盖所有发射分支，最终
  `pairVecDivVHCOuterIteration_nodes_size` 贯通 consume/emit/reinsert。
  因而 `nodes.size = divisor.size - 1` 可在全局 outer loop 中逐轮传递，
  无需把长度保持作为外加 oracle 假设。
- general heap division 的全局良基高次系数语义已闭合。
  `pairVecDivVHCOuterIteration_quotient_eq_of_frontier_lt_lead` 从真实 emit
  分支证明 frontier 低于 divisor lead 时单轮不追加 quotient；
  `pairVecDivVHCOuterLoop_quotient_eq_of_limit_le_lead` 沿 generated
  `degreeLimit` 良基递归证明其后所有实际循环均保持 quotient。
  `pairVecDivVHCOuterLoop_productAgreesAbove_lead_of_success` 随后对完整
  `OuterLoop` 作强归纳：lead 以上的 frontier 使用真实 consume/emit 系数桥
  并递归传递 canonical、ownership、coverage、node-size、homogeneous、heap
  order、cursor prefix、processed lead 与 consumed dividend；低于 lead 的
  frontier 仍执行到终止，但用 quotient 不变与 sparse-gap 定理闭合所有
  lead 以上系数。该定理准确刻画 quotient-with-remainder，而没有把任意
  非整除输入伪装成完整乘积等式。下一步在 SQF 调用点用真实整除前提消去
  低次余数并推出完整多项式等式。
- exact-division 余数消去与 raw selector 成功性已闭合。
  `sparsePolyZp_toPoly_degree_eq_head` 从 canonical sparse 首项的非零系数
  与严格降次链证明 L2 polynomial degree 等于数组首项次数；
  `pairVecDivVHCProduct_eq_of_agreesAbove_lead_of_dvd` 将全局高次系数一致
  改写为 `degree (dividend - quotient*divisor) < degree divisor`，再由 SQF
  调用点的真实 `divisor ∣ dividend` 推出该余数同时被 divisor 整除，故只能
  为零。这没有向执行器提供预期商。
  raw→safe 方向上，`pairVecDivVHCSelectFrontier_succeeds` 从非终止条件和
  heap chain ownership 推出 selector 必返回 `.ok`：非空 heap root 由 owner
  chain 保证是有效 active 节点，所有 dividend 读取则由分支边界保证。
  下一步继续闭合 consume/extract、emit activation 与 reinsertion 的成功性，
  最终组合为 OuterIteration/OuterLoop 的存在性定理。
- 单个 bucket-chain 节点的 raw→safe 成功边界已精确闭合。
  `pairVecDivVHCConsumeNode_succeeds` 从真实 node denotation 推出 node、
  quotient cursor 与 divisor cursor 的全部 checked bounds，并证明 advance
  与 exhausted 两条实际写回分支均返回 `.ok`。定理没有隐藏源码的关键
  assertion：当 cursor 恰好到达 quotient 末尾时，唯一剩余前提被明确收敛
  为 `nodeIndex = resetH`。因此后续 chain/root 成功性证明的真正任务已经
  精确定位为 exhausted rows 的顺序不变量，而不是泛化的“不会出错”假设。
- equal-degree consume 的 deferred-node 次数不变量已建立。
  新谓词 `PairVecDivVHCLinBelow` 记录本轮已推入 `lin` 的节点均已前进到
  严格低于当前 frontier 的乘积次数；空 `lin` 初始态直接成立。
  `pairVecDivVHCConsumeNode_preserves_linBelow` 跟随真实 advance/exhausted
  分支证明保持：advance 分支通过 node denotation 与 canonical quotient 的
  后继项严格降次得到新 monomial 低于 frontier，exhausted 分支则利用当前
  chain 节点不在旧 `lin` 中运输旧成员。该性质用于证明 `resetH` 行不可能
  已被本轮 deferred，从而可在 heap ownership 中找到它并排除越序耗尽。
- concrete VHC 行序的次数支配性质已闭合。
  `pairVecDivVHCNode_product_degree_gt_of_cursor_le_of_divisor_lt` 直接从
  quotient/divisor 两个 canonical 稀疏数组的严格降次顺序，以及两个 active
  节点各自的真实 denotation，证明：较早 divisor 行若 quotient cursor
  不更晚，其当前乘积次数必严格更高。该结论不依赖 heap 算法规格结果、
  预期商、fuel 或 L2 计算；它正是把源码 exhausted 分支的
  `nodeIndex == reset_h` assertion 从外加前提转化为可证明运行时不变量的
  数学核心。下一步将它与 `StateCovered`、`LinBelow` 和 selector heap-max
  组合，分别排除 `resetH` 行位于 lin 与 heap 的两种越序情形。
- selected frontier 上的耗尽顺序断言已从状态不变量推出。
  `pairVecDivVHCExhaustedNode_eq_resetH` 对任意当前 active frontier 节点
  证明：若其 quotient cursor 的后继恰为数组末尾，则节点索引必等于
  `resetH`。证明先由 `ResetReady` 排除节点位于已耗尽前缀，再由
  `StateCovered` 定位真实 `resetH` 行：若它在 deferred `lin`，行序次数
  支配与 `LinBelow` 矛盾；若它仍由 heap chain 拥有，则与 selector 的 heap
  最大次数矛盾。该结论完全消除了顺序 assertion 作为数学假设的必要性；
  下一步沿 ConsumeChain 的未访问集合运输入口快照，把此定理应用到每个
  尚未改写的真实 chain 节点。
- bucket chain 的 raw→safe 成功性已沿真实良基递归闭合。
  `pairVecDivVHCConsumeNode_preserves_base_reset_or_erased` 记录每步执行后
  entry `resetH` 要么尚未改变，要么对应节点已经从 `unvisited` 擦除；
  `pairVecDivVHCConsumeChain_succeeds` 以具体 `unvisited.card` 为终止度量，
  并运输“所有未访问节点仍等于 bucket 入口快照”。因此每次实际
  `ConsumeNode` 调用都能用入口状态的耗尽顺序定理消除 assertion，而不是
  使用步数预算或假定 raw 调用成功。`pairVecDivVHCConsumeRootBucket_succeeds`
  已将结论接到真实 root bucket 包装层。下一步组合 heap-size 良基的
  `ConsumeEqualDegree` 与 checked extract 成功性。
- generated heap extract 的 raw→safe totality 已闭合。
  `pairVecDivVHCSiftDown_succeeds` 从 concrete heap pointer validity 推出
  left/right/sentinel 三次 monomial 读取均成功，并沿执行器相同的
  `limit - child` 良基度量递归到实际选中的 child；每次覆盖 hole 后逐点
  运输 pointer validity。`pairVecDivVHCExtract_succeeds` 从非空 heap 的末项
  有效性实例化该定理，`pairVecDivVHCExtractChecked_succeeds` 再接通携带真实
  size-decrement 证明的 wrapper。整个证明不预设 extract trace 或成功结果；
  下一步可在 heap-size 良基的 `ConsumeEqualDegree` 中直接构造每轮 extract。
- deferred-node 的严格低次性质已沿完整 bucket chain 运输。
  `pairVecDivVHCConsumeNode_lin_subset_insert` 精确刻画单步后 `lin` 只可能
  新增当前节点；`pairVecDivVHCConsumeChain_preserves_linBelow` 同时沿 owner
  集合的擦除作良基递归，用 owner/lin 分离证明当前节点从未提前出现于
  deferred stack，并用 canonical quotient 的后继严格降次闭合 advance
  分支。`pairVecDivVHCConsumeRootBucket_preserves_linBelow` 已接通 root owner。
  因而 EqualDegree 的下一轮可继续使用同一 frontier 的耗尽顺序定理，而
  无需重置或弱化 `LinBelow`。
- selector 已可跨 root-bucket extract 精确运输。
  `pairVecDivVHCSelectFrontier_eq_of_root_degree_eq` 逐分支展开 generated
  selector，证明任何非空目标 heap 只要 root 次数等于原 frontier 次数，
  就会重新得到完全相同的 frontier record，包括 coefficient 与 dividend
  cursor。该桥覆盖 dividend 胜出、heap 胜出、次数 tie 与 dividend exhausted
  四类真实分支；EqualDegree 的 heap-size 递归因此可以在每个同次 root 上
  复用入口 selector 见证，而无需引入规格级 frontier。
