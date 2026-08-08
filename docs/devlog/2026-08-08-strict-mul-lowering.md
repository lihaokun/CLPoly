# HGCD 原始乘法 lowering

日期：2026-08-08

## 做了什么

- 固定 C++ `_classical_mul`、`_kar_mul`、`_mul` 三个方法的稳定 AST 哈希。
- 新增 schoolbook 乘法的 RawHeap 执行层。
- 内层严格使用源码的 `j_min..j_max` 点积区间，逐项读取 A/B，执行真实
  `_umul128` 和三字 `_add_carry3`。
- 外层对每个累加结果调用真实 `_lll_mod_preinv`，随后写入 C。
- 内外循环都以未处理区间为终止度量，并显式传播所有 RawHeap 故障。
- 从源码的 `j_min/j_max` 条件表达式证明每个 `A[j]` 与 `B[k-j]` 索引合法。
- 证明内层点积、外层逐系数写入及完整 schoolbook 入口在有效容量和非空输入
  下必然成功，并在所有写入后保持 RawHeap 分配布局。
- 定义与 C++ 点积循环读取完全相同 raw cells 的无界自然数求和执行，并证明
  它在有效索引下必然成功。
- 逐次组合已验证的 `_umul128/_add_carry3` 定理，证明最终三字累加器与
  “初始累加器 + 原始乘积和”在机器宽度 `2^192` 下同余。
- 由每个规范系数至多为 `p-1`，递归证明 raw 点积和不超过
  `项数 × (p-1)^2`。
- 使用项数小于 `2^64`、`p>1` 和既有 lazy-accumulation 容量不等式，证明
  点积和严格小于 `2^192`，从而把机器同余提升为普通自然数精确相等。
- 由精确点积容量推出三字高 limb 小于 p，直接应用已验证的生成
  `_lll_mod_preinv` 定理，证明实际归约输出等于 raw 点积和 `% p`。
- 进一步证明该输出转换到 `ZMod p` 后等于无界 raw 点积和的自然数转换。
- 定义与源码循环边界一致的 L2 系数递归和，并从 `SlicePolyRep` 的实际
  raw reads 逐项证明：无界 raw 点积和转换到 `ZMod p` 后，恰好等于同一
  区间上的 L2 系数乘积和。该桥没有调用 L2 乘法来代替 C++ 执行。
- 证明递归系数和等于 `j_min..j_max` 闭区间上的有限和；再用输入的
  `SlicePolyRep` 证明该区间外的每个标准卷积项均为零，最终得到源码点积
  区间精确等于 `Polynomial.coeff (left * right) k`。
- 对实际外层写循环证明帧性质：它只会改写 `C[k..lenC)`；任意不与这些
  cell 重合的 raw read 在完整剩余执行后保持不变。该性质将用于保持 A/B
  输入表示以及此前已经写出的 C 系数。
- 将单个源码系数的整条执行链闭合：实际 raw 点积、精确三 limb 累加、实际
  `_lll_mod_preinv`、`ZMod` 转换以及标准乘积系数全部相等。
- 从外层帧性质导出与 C 不同 allocation 的整个输入 prefix 保持不变，并
  证明 `CanonicalU64Prefix` 可沿相同 prefix 传递，为外层递归保持 A/B 的
  表示与规范剩余不变式。
- 建立外层循环的已完成输出前缀不变式：一次真实 `writeU64` 将已证明前缀
  精确扩展一个乘积系数，而剩余外层循环保持此前所有输出系数不变。
- 以 `lenC-k` 为良基度量闭合完整外层递归：每轮从当前 raw heap 执行真实
  点积和写入，扩展输出前缀，并通过 allocation 不相交性保持 A/B 的
  `SlicePolyRep` 与规范剩余；退出时 C 的全部 `lenC` 个 cell 均对应乘积
  系数。
- 从完整点态输出不变式构造实际 `SlicePolyRep`：先由 raw slice 的唯一观测
  表示逐系数识别缓冲区内内容，再用卷积 antidiagonal 与 A/B 声明长度证明
  `lenA+lenB-1` 之外的所有乘积系数为零。由此得到 `_classical_mul` 的完整
  raw slice 到 L2 乘积精化定理。
- 强化输出前缀不变式，使每个实际 `_lll_mod_preinv` 结果同时携带 `< p`
  的规范剩余证明，从而得到完整输出 `CanonicalU64Prefix`。
- 在素数模数下，从输入 `RawDensePolyRep` 的真实 normalization 结果推出
  A/B 末系数非零；用卷积最高项公式证明乘积末系数非零，再直接展开实际
  `normaliseU64` 的末 cell 分支，证明输出长度保持 `lenA+lenB-1`。
- 最终定理 `classicalMul_refines` 已将真实 `_classical_mul` 从两个规范 raw
  输入严格精化为 `RawDensePolyRep (left * right)`，包括终止执行、布局、
  规范剩余、完整 L2 表示和规范长度。
- 新增 `_kar_mul` 的完整 RawHeap 生成层：按源码顺序执行半区模加、奇数
  尾项复制、P0/P1/P2 三次递归、两轮原地模减、P0 `copyU64`、间隙清零和
  交叉项原地模加。递归以 `n` 为良基度量，并从 `n≥16` 证明 `n/2` 与
  `n-n/2` 都严格小于 `n`；没有引入执行次数参数。
- 对 Karatsuba 的三类原始循环分别建立终止执行与布局保持定理：半区相加
  循环验证四次实际读取和两次写入，交叉项模减循环验证读-读-写，最终组装
  循环验证偏移 C cell 与 P1 cell 的读-读-写；全部以剩余区间良基递归。
- 对三个 helper 建立精确帧定理：半区相加只可能改写 t1/t2 的剩余同索引
  cells，模减只可能改写目标剩余区间，组装只可能改写 `C[m+i]` 的剩余区间。
  每个定理逐次应用真实 `readU64_writeU64_ne`，可用于证明输入、递归结果和
  已完成输出在后续原地更新中保持不变。
- 对半区相加 helper 的任意当前索引证明最终 raw 值：真实四次 A/B 读取经
  `nmod_add` 写入 t1/t2 后，由帧定理保持到循环结束。进一步从输入
  `SlicePolyRep` 与 `CanonicalU64Prefix` 推出两个输出分别表示对应低/高半
  L2 系数之和，并且实际机器值继续严格小于 p。
- 将当前索引结论提升到完整剩余区间：每轮 t1/t2 写入后，利用 scratch 与
  A/B allocation 不相交性逐 cell 构造 `SameU64Prefix`，重新传递输入
  `SlicePolyRep` 和规范剩余，再对 `m-(i+1)` 做良基递归。由此从 i=0 可
  得到全部 `k<m` 的 t1/t2 raw 值与 L2 半区和语义。
- 定义有限 L2 半区和多项式并证明其系数公式；再以 t1/t2 全部 raw cell 的
  唯一观测把区间结论提升为两个完整 `SlicePolyRep`，同时从每个实际
  `nmod_add` 输出 `<p` 得到两个 scratch slice 的 `CanonicalU64Prefix`。
  该定理覆盖长度 `2m` 的公共部分，奇数 n 的第 h 个尾项仍按源码分支另证。
- 将源码内联的 `if (h>m)` 尾项块原样抽为 RawHeap helper `karOddTail`，主
  `_kar_mul` 保持相同调用位置与故障传播。证明该块在 A/B 长度 `m+h`、
  t1/t2 长度 h 下终止并保持布局；同时证明其两次写入保持所有与 t1/t2
  allocation 不相交的 raw prefix，为输入和既有 scratch 前缀传递提供帧。
- 在奇数分支逐条识别尾项执行：证明 A/B 的 `2m` 实际 reads 经 t1/t2 第 m
  个 cell 的两次原样 writes 后保持到最终 heap；再由完整输入
  `SlicePolyRep` 和规范剩余得到这两个 raw 值分别等于原 L2 的第 `2m`
  个系数且严格小于 p。
- 证明尾项块保持 t1/t2 自身的 `[0,m)` raw prefix；将已有长度 m 的半区和
  `SlicePolyRep` 传递到最终 heap 后，使用真实第 m 个尾 cell 分别扩展为
  长度 `h=m+1` 的完整 slice。输出多项式精确为公共半区和加第 `2m` 系数
  的单项式，并逐分支证明整个新 slice 仍为规范余数。
- 新增通用 raw prefix 提取桥：从较长 `SlicePolyRep` 的相同底层 allocation
  构造任意较短 prefix 的唯一 L2 表示，并证明前缀内系数与完整输入一致。
  该桥用于奇数 n 时从 `m+h` 输入安全取得公共 `2m` 部分。
- 将公共半区相加与奇数尾项按源码顺序组合为 `karPrepareHalves`，主
  `_kar_mul` 直接调用该真实 RawHeap 块；证明在 `m≤h` 和完整容量下该组合
  必然成功并保持布局。
- 证明完整准备块保持所有与 t1/t2 allocation 不相交的 raw prefix：先用
  公共循环的逐写帧传递到中间 heap，再用尾项块的区域帧传递到最终 heap。
  因而完整 A/B 表示可安全跨准备阶段保持。
- 证明 `karHalfSumPoly` 只依赖输入前 `2m` 个系数：任意与完整多项式在该
  prefix 上逐系数一致的截断表示，计算出的半区和多项式与完整输入版本相等。
- 闭合统一准备精化 `karPrepareHalves_refines`：从完整 `m+h` 输入提取真实
  `2m` prefix，执行公共半区相加，保持完整 A/B 表示到中间 heap，再按
  `h=m` 的无写尾分支或 `h=m+1` 的真实尾项分支组合。最终 t1/t2 都获得
  长度 h 的 `SlicePolyRep (karPreparedPoly ...)` 与完整规范剩余证明。
- 定义源码递归真正可达的精确 scratch 需求 `karScratchNeed`：基例为 0；
  递归层为 t1/t2/P0/P1 当前分区之和，加共享 `recScratch` 中两个子规模需求
  的最大值。三个乘法调用顺序共享同一区域，因此这里不能错误地相加。
- 证明源码二分恒为 `h=m` 或 `h=m+1`，且在递归分支 `n≥16` 时 m/h 都
  严格小于 n。后续容量证明将使用精确递归需求，避免奇数分割下并不成立的
  “rec 尾部单独拥有 6h cells”假设。
- 从精确递推式证明当前 t1/t2/P0/P1 四个分区总长不超过
  `karScratchNeed n`，并证明共享递归尾区的 `max` 同时支配 m/h 两个子调用
  的精确需求。全局 `karScratchNeed n≤6n` 尚需携带递归余量的不变式；未用
  奇数层不成立的朴素 `6h` 归纳替代。
- 进一步审计发现 `6n` 不只是难证，而是在连续奇数分割下为假：例如
  `n=2^31+1` 时精确需求已超过 `6n`。因此修正 C++ 容量合同：公开 mul 的
  Karatsuba scratch 从 `6n` 增至 `7n`，raw `_mul` 总 scratch 从 `7n`
  增至 `8n`，HGCD 工作区同步调整；设计文档和 `_gcd_hgcd` AST 锁已更新。
- 以 n 为良基度量证明 `karScratchNeed n ≤ 7n`。递归步分别使用 m/h 子结论、
  共享区的 max 上界以及 `h=m ∨ h=m+1`，覆盖全 Nat 范围。
- 从一块长度为 `karScratchNeed n` 的真实 allocation 推导源码指针算术切出的
  t1、t2、P0、P1 和共享 recScratch 五段全部有效；嵌套 `ptr.add` 的地址等式
  也在 RawPtr 模型中逐项证明，没有把容量不等式当作指针有效性的替代。
- 新增 UInt64 slice 的地址级不相交谓词及对称、缩短、不同 region、相邻区间
  和同基址偏移区间引理。旧证明只用 region 不等表达不别名，无法描述 C++
  同一 scratch allocation 内的 t1/t2/P0/P1；后续 helper 与 schoolbook 帧定理
  将改用这一谓词，避免为套用旧定理而虚构独立 allocations。
- 将 Karatsuba 半区准备的整条语义链改为地址级 slice 不相交：公共相加、奇数
  尾项、自身前缀保持、完整 slice 构造及统一准备定理现在都直接接受同一
  scratch allocation 内相邻的 t1/t2，不再要求事实上不成立的 region 不同。
- schoolbook 的 coefficient-prefix 与 slice 精化同样推广到地址区间不相交；
  原公开 raw 定理仍可从不同 allocation 推出该条件，而递归 P1 基例可直接
  使用 scratch 内的偏移区间。
- 闭合 `karMul_ok`：按生成 L1 的真实控制流依次执行半区准备、P0/P1/高半区
  三次递归乘法、两次模减、P0 copy、边界补零和交叉项组装。每次调用都从
  精确 scratch 分区与 SameLayout 传递实际 slice 有效性；基例执行已证明的
  schoolbook。递归只以 n 良基下降到 m/h，没有 fuel 或 L2 回退。
- 证明 Karatsuba `karSubLoop` 与加减生成层的 `subCommonLoop dst dst sub`
  在每次 read/read/write 及错误传播上结构相同；该等式只复用 raw 执行结构，
  不调用带规范化的高层 `_poly_sub`，也不替代 Karatsuba 的实际 heap 执行。
- 直接在地址级不相交条件下证明 `karSubLoop_value`：原目标 cell 与 sub cell
  穿过更早 dst 写入保持不变，当前 cell 写入真实 `nmod_sub`，再穿过其余循环
  保持到最终 heap。由逐 cell 结论闭合完整 `SlicePolyRep (left-right)` 和输出
  `CanonicalU64Prefix`，可用于同一 scratch allocation 内相邻 P1/P0 分区。
- 闭合 `karAssembleLoop_value`：更早的 C 写入保持目标 `C[m+k]` 与 P1[k]
  原始 reads，当前迭代执行真实 `nmod_add` write，后续写入再保持该结果。C/P1
  的帧条件采用地址级 slice 不相交，可直接对应输出 allocation 与 scratch。
- 将组装逐 cell 结论提升为完整多项式表示：低于 m 的 C prefix 由真实帧保持，
  `[m,m+count)` 逐项加入 P1，slice 外双方系数均为零，最终精确得到
  `basePoly + X^m * crossPoly`；同时证明整个 C 输出区间仍为规范模 p cells。
- 新增 raw slice 二分桥：从长度 `lowLength+highLength` 的同一表示构造底部
  prefix 与 `ptr.add lowLength` 尾部的唯一 L2 多项式，并逐系数证明原多项式
  精确分解为 `low + X^lowLength * high`；尾段系数通过真实 `readU64_add`
  地址等式关联，未把数组切片作为未证假设。
- 证明奇偶二分下 `karPreparedPoly` 精确等于 `low+high`。偶数分支覆盖 m 个
  公共和，奇数分支额外证明源码尾 cell 对应 high 的第 m 项；所有 slice 外
  系数均由表示长度定理归零。
- 在 Polynomial 环中闭合 Karatsuba 恒等式：P0、`(low+high)` 乘积减 P0/P2
  所得交叉项、以及移位 `X^(2m)P2` 的和精确等于原左右输入之积。该恒等式
  将作为递归 raw 输出组合到最终 C 表示的代数终点。
- 将 `karOddTail` 与完整 `karPrepareHalves` 的 guard 帧定理从 region 不同推广
  为任意地址级不相交 slices。每次 t1/t2 写入都用实际写索引和 guard 读索引
  证明地址不同，因此同一 scratch allocation 内的 P0/P1/recScratch 可作为
  guard 被可靠保持。
- 新增 schoolbook 完整调用的地址级帧定理：从生成 `_classical_mul` 的成功
  执行还原真实 outer loop，并逐 guard cell 应用已有逐写帧，得到任意与 C
  输出 slice 不相交区域的 `SameU64Prefix`。这是 Karatsuba 递归基例的帧桥。
- 证明总 slice 与 guard 地址不相交可下推到任意合法 `ptr.add start` 子 slice；
  该引理把顶层 C/scratch 帧合同系统地传给 t1/t2/P0/P1/recScratch，而无需
  为每一对指针重复展开地址算术。
- 建立 `SameU64Prefix` 的跨 heap 传递，以及单次 raw write、真实 copyU64、
  完整 karSubLoop、完整 karAssembleLoop 的 slice 级帧包装。每个包装仍调用
  已证明的逐 read/write 定理；它们将用于按生成控制流组合递归 frame。
- 证明 UInt64 指针嵌套偏移满足 `(ptr.add a).add b = ptr.add (a+b)`，并将
  顶层 `scratch[0..karScratchNeed n)` 与任意 guard 的不相交性一次性下推为
  t1、t2、P0、P1、共享 recScratch 五段各自的不相交合同。证明使用精确需求
  递推式验证每段范围，并显式归一化源码的嵌套 pointer arithmetic。
- 强化 `karMul_ok` 的同一个 n 良基归纳，使每层除成功执行与 SameLayout 外，
  还返回任意同时避开输出 C 和该层 scratch 的 guard `SameU64Prefix`。基例
  使用真实 schoolbook frame；递归层依次组合半区准备、P0、P1、高半区、
  两次模减、copy、补零与组装九个实际 heap 阶段。
- 三个递归调用的 frame 直接取自各自良基归纳结果，recScratch 的 child 需求
  由 max 缩短，输出段由五分区不相交定理取得；后处理分别使用已证明的
  sub/copy/write/assemble slice frame。新增 `karMul_preserves_prefix` 将任意给定
  成功运行与唯一 ok heap 对齐后导出该完整帧性质。
- 消除半区准备语义链剩余的四个 region 假设：t1/t2 与 A/B 的关系全部改为
  地址级 slice 不相交。公共相加递归在每次 write 后以实际 i/j 地址保持 A/B，
  完整准备层再把 h 长度合同缩短到 m、把输入长度缩短到 2m 后传入。
  因此 P1 子递归中 t1/t2、其 recScratch 与上一层 scratch 同属一个 allocation
  时仍可真实应用准备精化，不需要虚构不同 region。
- 强化 schoolbook slice 精化，使其除乘积 `SlicePolyRep` 外显式返回整个输出
  区间的 `CanonicalU64Prefix`；该性质直接来自每个真实点积模约减写入形成的
  `ClassicalCoeffPrefix`，递归基例输出可安全进入后续 nmod 运算。
- 新增真实 `copyU64` 的 slice+canonical 联合精化：利用 copy 的逐 cell 内容
  定理把源值与目标 read 对齐，同时传递源 `<modulus` 条件。P0 copy 到 C 后
  因此保留完整 L2 表示和机器规范性，而非只证明 memcpy 成功。
- 建立相邻 raw slices 的反向拼接定理：长度 low/high 的表示可在同一连续
  allocation 中合成为 `low + X^lowLength * high`；证明通过完整 slice 的唯一
  表示与已证 split 桥对齐。相应规范前缀逐索引使用 `readU64_add` 拼接。
- 闭合末尾补零写入：真实 `writeU64 ptr length 0` 保持旧 prefix，新增 cell
  读取恰为零，并把 `SlicePolyRep` 从 length 扩展到 length+1 而多项式不变；
  同时由非零 modulus 证明新零 cell 仍规范。这覆盖 P0 后的源码分隔零。
- 闭合 Karatsuba 后处理的 copy/zero 基段：从真实 `copyU64 C P0
  (2*m-1)` 和随后真实零写入出发，同时保持已在 `C+2*m` 中的 P2
  高位 slice，最终拼接出 `P0 + X^(2*m)*P2` 的完整 `SlicePolyRep`
  和规范前缀。整段证明逐次使用 copy/write 的实际 heap 与地址帧，
  没有用一次性规格替换生成控制流。
- 将真实递归 `copyU64` 的内容语义从“不同 allocation region”推广到
  逐地址不相交。证明在每次头 cell 写入后保持源尾段，然后对同一
  allocation 的相邻尾 slices 继续递归，并把逐 cell 内容提升到
  `SlicePolyRep` 与 canonical 结果。因此子层 Karatsuba 的 C/P0 都在父层
  scratch allocation 中时也可直接应用，不再需要虚构 region 不等。
- 证明每层 Karatsuba scratch 中 t1、t2、P0、P1、recScratch 五个
  连续区间的成对地址不相交，并将规范 offset 形式对齐到生成 L1
  中字面出现的嵌套 `ptr.add`。这些合同同时支撑 P0/P1 递归输出与
  共享子 scratch、两次交叉项减法以及 P0→C copy，不依赖 allocation
  region 不等。
- 将真实 `karPrepareHalves` 提升为递归子调用可直接消费的接口：
  除了返回 t1/t2 对 `karPreparedPoly` 的 slice 表示和 canonical 性，还返回
  真实执行的 `SameLayout`，并由每次 t1/t2 write 的地址帧证明 A/B 的
  `SlicePolyRep` 与 canonical 性在结果 heap 中保持。因此紧随其后的 P0
  递归不需要假设输入数组没有被改写。

## 当前边界

本步只完成 schoolbook 执行层。Karatsuba 与 `_mul` 的源码身份已经锁定，
但其零填充、scratch 布局、三次递归乘法、交叉项减法和组装尚未 lowering，
因此不能称为 `_mul` 精化完成。下一步先证明 schoolbook 的终止执行、点积
内容和完整输出表示已闭合。Karatsuba 的 raw 生成执行和良基递归现已落下；
下一步证明各 helper 的终止/帧性质、scratch 分区不相交性及 Karatsuba 恒等式。
在这些语义精化和 `_mul` 分派完成前仍不能称为 `_mul` 精化完成。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictMul.lean`
- `proof/lean/CLPoly/Impl/StrictMulRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_mul_source.py`
- `docs/devlog/2026-08-08-strict-mul-lowering.md`
