# HGCD raw lowering 起点

日期：2026-08-08

## 本步完成

- 确认通用 cpp2lean 对 `_hgcd_iter` 的当前输出仍含二级指针解引用的
  `sorry`、`partial def`、未建模 `swap` 与 `sizeof`，该输出没有进入严格边界。
- 从 HGCD 真实调用链底部 `_mat_one` 开始建立专用 RawHeap lowering：依次执行
  C++ 的 `M.poly[0][0] = 1`、`M.poly[3][0] = 1`，并更新四个长度字段。
- 建立 bounds-safe 的四项矩阵访问器与 `HgcdMatPolyRep`，证明真实写入后的
  矩阵精确表示 `[[1,0],[0,1]]`，同时保持 heap allocation layout。
- `_mat_one` AST 哈希已加入 HGCD source gate；gate 同时扫描生成层、raw
  精化层和代数层的禁止构造。
- 将 `_mat_row_update` 的二级指针引用消除为显式返回状态，但保留真实
  `_mul`、`_poly_add`、乘积长度计算以及三个指针/长度交换。先闭合源码的
  `Q == 0 || M[i0] == 0` 分支：证明它不访问 heap，只交换两项矩阵描述符，
  并保持 `HgcdMat.Valid`。该函数的 AST 哈希和 `_mul`/`_poly_add` 调用结构
  也已加入 source gate。
- 补齐非零分支的成功执行分解：若生成函数返回成功，则证明其中同一参数的
  `_mul` 必先成功，随后同一 heap 上的 `_poly_add` 必成功，并精确恢复
  product length、返回 T/t 指针和更新后的四项矩阵 descriptor。该引理排除
  后续证明绕开这两次真实调用、只按规格构造矩阵结果的可能。
- 强化成功分解以证明 descriptor 的精确落点：在 `i0 != i1` 下，新 `i0`
  指向加法输出 t/`sumLen`，新 `i1` 指向旧 `i0`/旧长度。随后将真实
  `_poly_add` 的 `RawDensePolyRep` 定理接到该落点，得到新行项精确表示
  `old[i1] + product`。其中 product 将由严格 `_mul` 定理实例化为
  `Q * old[i0]`，而不是在矩阵层重新计算。
- 已将该 product 前提闭合到真实 dispatcher：按 C++ 的长度条件选择参数
  顺序，覆盖 schoolbook 与良基 Karatsuba，并用执行确定性确认 dispatcher
  返回的就是行更新所观察到的同一个 heap。交换参数后仅使用乘法交换律统一
  为 `Q * old[i0]`；没有另造规格执行或 L2 回退。
- 同一真实 `_mul` 的写区域 frame 现已证明旧 `i1` 项逐 cell 不变；由此
  重新导出其有效切片、canonical 系数、L2 多项式表示和 normalization，
  不再把“乘法后旧项仍可供 `_poly_add` 使用”留作语义假设。
- 严格 `_poly_add` 现有完整外部 allocation frame；HGCD 用它证明加法后
  旧 `i0` 缓冲区仍保持 normalized polynomial 表示，因此 descriptor 交换
  到新 `i1` 后的内容语义也不再缺失。
- 非零 `_mat_row_update` 已组合成单一端到端定理：从同一个生成函数成功
  执行推出 `new[i0] = old[i1] + Q * old[i0]` 且
  `new[i1] = old[i0]`，两项均为最终 heap 上的 `RawDensePolyRep`。

## 下一步

`_hgcd_iter` 已建立完整 reference-eliminated RawHeap lowering：保持 C++ 的
`_mat_one`、A/B 两次 memcpy、真实 `_poly_divrem`、指针旋转和两次
`_mat_row_update`。循环直接按当前 `lenB` 良基递归；下降证据由实际 divrem
返回的 `lenR < lenB` 控制流定理提供，不使用 fuel 或额外运行时断言。

初始化前缀现已抽成 `hgcdIterInit`，完整入口严格等于初始化成功后进入同一
良基循环；真实 recursive memcpy 也已有完整 `RawDensePolyRep` 传递定理，
包括有效性、canonical、L2 slice 语义和 normalization。
- `_mat_one` 的两次真实 `writeU64` 现已具备外部 prefix frame，因此后续
  初始化组合可证明单位矩阵写入不会改变 `a/b` 输入，而不是假设输入幸存。
- `_mat_one` 精化结果现还显式给出 `poly` 指针数组不变且 `len` 精确为
  `[1,0,0,1]`，使后续 memcpy 的矩阵非别名合约可直接由原始 allocation
  条件实例化。
- 四项 `HgcdMatPolyRep` 现可由逐项 layout/prefix frame 整体传递，并已
  专门连接到真实 recursive memcpy；矩阵切片有效性保持为显式 L1 条件。

下一步证明初始化 copy 的表示传递，以及每轮 divrem/两次行更新共同保持
Euclid 对与矩阵变换关系，最终导出循环停止条件和 GCD 不变式。

`_hgcd_recursive` 的早退矩阵复制现已建立长度描述符不变式：循环入口 `i`
之前的槽位已与递归结果 `R` 相同，每次真实 `copyU64` 后同步执行对应的
`Array.set`，成功完成四项后返回矩阵的整个 `len` 数组精确等于 `R.len`。
递归度量仍为 `4 - i`，没有 fuel、规格执行或 L2 回退。该定理已加入严格
HGCD source gate；下一步是将同一循环的四次 heap 写入组合成每个矩阵项的
最终 `RawDensePolyRep`，从而闭合早退矩阵的完整语义。

早退矩阵的四项 heap 内容语义现已闭合。新增的 refine workspace 明确要求
每个目标与全部源项不相交、四个目标两两不相交；归纳证明逐次调用真实
`copyU64` 的 `RawDensePolyRep` 精化定理，并 frame 所有源项和此前已复制的
目标项。入口定理同时给出返回 `len = R.len` 以及四个返回矩阵项精确表示
`R` 的四个 L2 多项式。该证明没有从描述符相等推测 heap 内容，也没有调用
任何规格 oracle 或 L2 执行回退。

矩阵复制循环另有通用 raw frame 定理：只要一个活跃多项式切片与四个目标
区域逐项不相交，真实循环成功后其完整 normalized 表示保持不变。该结论将
直接用于外围早退函数的 A/B 两次恢复，保证随后的可选矩阵复制不会覆盖
刚写回的 A、B。

`hgcdRecursiveEarlyReturn_refines` 现已端到端组合早退块：先由真实
`copyU64 A a2` 得到重构后的 A，再 frame `b2` 和递归矩阵源；随后由真实
`copyU64 B b2` 得到 B，并 frame A 与矩阵源。`computeM=true` 时继续调用
同一个四项矩阵循环及其内容精化，同时利用通用 frame 定理保持 A/B；为
false 时证明返回矩阵就是原 M。返回的 `lenA/lenB/sgn` 也精确对应 C++
赋值。所有别名条件集中在纯物理 workspace 中，没有把预期 L2 值藏进执行
合约。

非早退路径现已继续降低到中间 divrem：`hgcdRecursiveMiddle` 调用真实
`dense_upoly_zp__poly_divrem_ir(q,d,a2,b2,...)`，随后严格按源码计算
`k = 2*m-lenb2+1`、`c0=b2+k`、`d0=d+k` 及两个截断长度。精化定理从
同一次 raw 执行导出商余式恒等式、余式次数下降、`lenD < lenB2`、商余式
canonical/normalized 表示和第二次 HGCD 输入布局。这一步没有把 divrem
替换为 EuclideanDomain 的规格计算。

第二次 HGCD 调用的良基算术也已闭合：在重构输出长度不超过外层 `lenA`
的物理界下，源码 `k` 使 `lenC0 < lenA`；若 `c0` 非空，则同一偏移作用于
divisor/remainder，并由真实 divrem 的 `lenD < lenB2` 推出
`lenD0 < lenC0`。中间精化定理现直接返回第一条下降事实，供最终
`_hgcd_recursive` 的 `termination_by lenA` 使用。

middle 精化输出已进一步强化为最终 heap 上的三个完整
`RawDensePolyRep`：商 `q`、余式 `d`，以及被真实 divrem 保持不变的除数
`b2`。因此第二次调用的 `c0=b2+k`、`d0=d+k` 将从真实 raw 表示切分，
不需要按商余式 L2 等式重新制造输入缓冲区。

第二次调用的 raw 输入切分现已闭合：在源码确实进入非空 `c0` 分支时，
由 `lenC0` 的生成公式反推出 `k ≤ lenB2`，直接对最终 heap 中保持的 `b2`
执行表示切分；`d0` 在 `k ≤ lenD` 时同样从真实余式切分，否则由工作区
提供的零长度有效性得到零多项式表示。整个桥不执行写 heap，也不调用 L2
算法生成 `c0/d0`。

为无额外运行时 guard 地组合递归内 iterator 分支，生成层现已证明完整的
descriptor 有效性链：`_mat_one`、每次 `_mat_row_update`、良基
`hgcdIterLoop`、初始化以及完整 `_hgcd_iter` 的每个成功结果都保持
`HgcdMat.Valid`。循环证明与执行使用相同的 `state.lenB` 终止度量，下降仍
直接来自真实 divrem 返回的余式长度。

递归函数的 iterator arm 现已组合为单一生成执行：先调用真实 `_hgcd_iter`，
再以其已证有效的矩阵执行两阶段稳定化，最后执行带交叉保护的 `pA/pB`
写回。成功分解定理强制暴露这三次同序执行及所有返回字段，防止后续精化
绕开稳定化或 memcpy。该组合没有增加源码外的动态 descriptor guard。

初始化 copy 的组合现已闭合为 `hgcdIterInit_refines`：从一次真实初始化执行
同时导出 identity matrix、`A = a`、`B = b` 的最终 heap 表示，以及 A/B/T/t、
长度和初始符号的精确状态。所有 copy/矩阵/input 非别名与切片有效性均为
显式 L1 前置条件。

递归 iterator 输出归一化的矩阵 frame 也已闭合：无论源码进入交叉别名保护
分支，还是分别跳过/执行 `a3`、`b3` 的条件复制，每次真实 `copyU64` 都保持
四个稳定矩阵项的 raw 表示，且组合后给出初末 heap 的 `SameLayout`。该证明
只依赖显式切片有效性与目标/矩阵项不相交，不对 heap 或矩阵内容作规格替换。

矩阵稳定化现在也具有完整的 live-polynomial frame：证明分别沿生成的 staging
与 restore 良基四步循环展开，每次由真实 `copyU64`、目标有效性和区域分离
保持任意 `RawDensePolyRep`，并组合得到初末 heap 的 `SameLayout`。因此递归
iterator 产生的 A/B 可穿过稳定化进入最终输出复制，无需把其 L2 内容作为
workspace 字段或重新计算。

递归 iterator arm 已端到端闭合：同一个成功的生成执行被强制分解为真实
`_hgcd_iter`、真实矩阵稳定化和真实别名敏感输出复制；组合定理返回最终
`a3/b3` raw 表示、稳定矩阵四项 raw 表示、输入输出线性变换、与源码 `sgn`
同步的行列式、GCD 不变式以及 `lenB < lenInputA/2+1` 停止界。新增 finalize
workspace/provider 只含容量、指针相等桥和分离条件，不含任何 L2 多项式值。

两处 A/B 重构共用的完整源码顺序现已生成化：`ReconstructB`、B 的
`LiftHigh`、`ReconstructA`、A 的 `LiftHigh` 被组合成一个 raw helper；成功
分解定理固定这四次执行的顺序，并将返回长度精确绑定到两次真实
`normaliseU64` 的结果。以 `shift=m,R` 实例化第一处，以 `shift=k,S` 实例化
第二处，不引入新的算法或规格执行。

最终矩阵乘法依赖已开始严格降低。`hgcdMatMulEntry` 依次执行两个真实 guarded
乘积，并逐字保留源码的四个尾分支：双非零 `_poly_add`、仅 PQ、仅 RS 的
`copyU64`、双零。`hgcdMatMulLoop` 再以 `4-i` 为度量按行列索引公式执行四
项，逐项安装实际返回长度并证明最终 descriptor 有效；没有直接构造 L2
矩阵乘积。

`hgcdMatMulEntry_refines` 已将单项 raw 执行闭合为 `P*Q + R*S`。两个乘积
均来自真实 guarded `_mul`；第一项跨第二次乘法由 prefix frame 保持。四个
尾分支分别用真实 `_poly_add`、保持 PQ、真实 `copyU64` 或零长度表示完成，
零乘积都从实际返回长度为零及 `SlicePolyRep` 推出。物理 provider 仅含容量、
别名和分离条件。

HGCD 迭代的矩阵长度不变式现已贯穿完整的良基循环。证明先顺序消费真实
`_poly_divrem`、row `(2,3)` 更新和 row `(0,1)` 更新，取得商长度界、余数
严格下降以及四个 descriptor 长度界；随后用这些界证明下一状态仍满足四项
互补长度不变式，并以源码的 `state.lenB` 严格下降递归。迭代器的初始化、
停止分支以及递归 iterator arm 的 stabilization/store 也已接通，因此返回
矩阵四项相对于原输入长度的界和 `result.lenB ≤ result.lenA` 都来自实际执行。

迭代循环现还同步保持 `final.lenA ≤ inputLength` 和 `0 < final.lenA`；前者由
每步 `next.lenA = old.lenB ≤ old.lenA`，后者由实际循环 guard 保证递归步的
旧 B 非空。结合返回矩阵第 0/2 项的互补长度界，新的算术闭合定理已证明第一
次成对重构的三个真实长度分量——移位高部、`R[2]*a_lo`、`R[0]*b_lo`——
均不超过原 `len_a`，从而得到重构结果 `lenb2 ≤ len_a`。这正是中间真实
`_poly_divrem` 已有严格余数下降定理的缺失前提。

中间真实 divrem 现进一步直接返回第二次 HGCD 调用所需的全部长度事实。
提前返回 guard 失败且 `m>0` 时，源码公式 `k=2*m-lenb2+1` 保证 `lenc0>0`；
同一次 divrem 的 `lend<lenb2` 在共同切去 `k` 后给出 `lend0<lenc0`，而第一
重构界使已有定理给出 `lenc0<len_a`。因此第二次调用既满足输入长度严格有序，
也对外层 `len_a` 良基度量严格下降，不需要运行时计数器或额外递归假设。

完整 `_hgcd_recursive` 源码函数体已经组合为 `hgcdRecursiveBody`。它固定
`HGCD_CUTOFF=100`，按源码顺序连接 base、第一次 iterator/递归分派、第一次
成对重构、提前返回、中间 divrem、第二次 iterator/递归分派、最终成对重构和
可选矩阵合并。两处递归调用目前显式参数化为同一个调用接口，除此之外没有
未执行的步骤；各分支先通过真实执行定理取得矩阵有效性，再转换为统一返回
记录。下一步用两处已证明的 `lenA` 严格下降把该接口封装为良基不动点。

第一次递归切片的其余良基前提也已直接从源码指针算术闭合：共同切去
`m=lenA/2` 后，原始 `lenB<lenA` 保持为 `lenB0<lenA0`，且 `lenA0>0`；
配合非 base guard，已有 `lenA0<lenA`。因此第一次递归调用现在具备完整的
正长度、严格有序和外层度量下降三项证明。

最终矩阵合并块也已按源码顺序闭合：`hgcdRecursiveCombineMatrix` 先实际执行
两列商矩阵更新，再把其返回的 matrix/heap 直接交给完整四项 `_mat_mul`。
对应 raw 定理从两次生成调用分别取得 `[[q,1],[1,0]]*S` 与
`R*([[q,1],[1,0]]*S)`，没有用单独的 L2 赋值替代执行。对左矩阵的跨调用
保持通过物理 layout/prefix frame 表达，provider 不携带预期多项式结果。

第二次高半部 HGCD 返回后的完整尾部现由 `hgcdRecursiveFinish` 固定：它先按
源码顺序执行 B/A 成对重构，再仅在 `compute_M` 为真时把同一返回 heap 交给
最终矩阵组合，最后原样返回 `-(sgnR*sgnS)`。执行分解定理分别钉住
`compute_M` 两个分支，因而良基主函数不能跳过重构或用规格矩阵替代可选调用。

成对重构的 raw 语义现已闭合。`hgcdRecursiveReconstructPair_refines` 逐一消费
真实的 `reconstructB`、B 的零填充/加高半部/归一化、`reconstructA`、A 的
零填充/加高半部/归一化，并返回两个精确的 L2 多项式。B 跨后两次写操作的
保持以及矩阵项/低半部/高半部跨调用的保持均由 heap layout/prefix frame
给出；workspace/provider 不含最终 A、B 或任何预期多项式。

重构长度也不再只使用较松的 `2*max-1` 缓冲容量。生成层新增从真实 guard
直接推出的 `mulTerm.length ≤ lenLeft+lenRight`；单次 A/B 重构再通过实际
`_poly_sub` 的归一化长度界得到两个乘积长度之最大值，成对重构最后通过
`liftHigh` 的真实归一化界得到
`lenB ≤ max(shift+highLen, max(r2+lowA, r0+lowB))`。因此主递归只需从
HGCD 不变式证明三个分量均不超过原输入长度，即可取得 `lenB2 ≤ lenA`，
无需加入运行时 `hdec`。

矩阵 descriptor 的递归长度不变式现已建立初始端。生成层从真实 `_mat_one`
和两次 copy 的返回值证明初态矩阵长度恰为 `[1,0,0,1]`、操作数长度仍为
输入长度。raw 层据此建立四个互补界：每一行中与当前 A/B 配对的矩阵项
长度之和不超过 `inputLength+1`。下一步沿真实 divrem 与两次 row update
保持该不变式，以导出第一次 HGCD 返回的 `R[0]`、`R[2]` 分量界。

单次 row update 的 descriptor 长度界已从真实执行闭合。零分支精确读取源码
的两项交换；非零分支从同一次 `_mul` heap 过渡取得 product 的有效切片，再
对同一次 `_poly_add` 使用归一化长度界。因此更新项长度不超过
`max(oldOther, lenQ+oldCurrent-1)`，交换项等于旧 current，其余两项保持；
没有把 sum 长度或 descriptor 结果作为 workspace 假设。

完整四项 `_mat_mul` 的 raw 语义现已闭合。新的单项 frame 定理按真实的
两次 guarded 乘法和 add/copy/skip 尾分支，保持两个输入矩阵以及先前已完成
的输出项。`hgcdMatMulLoop_refines` 再沿生成的 `4-i` 递归逐项执行
`hgcdMatMulEntry`，最终返四个精确的 L2 矩阵乘法项。workspace/provider
只保存容量和分离事实，不含预期多项式或矩阵乘积。

第二次递归后的 `S_modified = [[q,1],[1,0]] * S` 源块已开始独立降低。
`hgcdMatSwapRows` 精确表示四次 descriptor swap；随后两次
`hgcdMatQuotientEntry` 依次对两列执行源码中的非零 guard、按长度选择的真实
`_mul`、原地 `_poly_add` 及单个长度更新。`hgcdMatApplyQuotient_exec` 固定两次
实际执行的顺序和中间 heap/matrix，不复用会产生不同 descriptor swap 的
`_mat_row_update`。下一步是为这两次执行建立完整 raw 矩阵语义与 frame。

商矩阵更新的 raw 语义现已闭合。单列定理在非活跃分支从真实零长度表示
推出乘积为零；活跃分支则从同一次生成 `_mul` 和原地 `_poly_add` 推出
`top + q*bottom`，同时逐项 frame 另外三个矩阵项和商 `q`。完整定理先
证明四次 descriptor swap 只重排 raw 表示，再顺序消费两次真实列更新，
得到 `[[q,1],[1,0]] * S` 的四项 raw 语义。两层 workspace/provider 均只含容量、
别名和分离条件。
