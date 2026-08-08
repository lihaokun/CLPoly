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

初始化 copy 的组合现已闭合为 `hgcdIterInit_refines`：从一次真实初始化执行
同时导出 identity matrix、`A = a`、`B = b` 的最终 heap 表示，以及 A/B/T/t、
长度和初始符号的精确状态。所有 copy/矩阵/input 非别名与切片有效性均为
显式 L1 前置条件。
