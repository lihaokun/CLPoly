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
