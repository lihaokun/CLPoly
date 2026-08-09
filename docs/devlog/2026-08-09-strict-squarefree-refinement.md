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
