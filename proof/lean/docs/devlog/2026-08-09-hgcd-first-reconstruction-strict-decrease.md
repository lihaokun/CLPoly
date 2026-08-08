# HGCD 首次重构严格下降

## 结论

首次递归 HGCD 返回后，真实的成对重构结果满足
`reconstructed.lenB < lenA`，因此可直接作为后续良基递归边的下降证据，
不需要 fuel、运行时下降保护或 L2 回退。

## 证明来源

- C++ 调用入口已有的严格次序 `lenB < lenA`；
- 首次 HGCD 返回值的停止界与 `lenA > 0`；
- 四个真实矩阵描述符中的第 0、2 项长度界；
- `_mul` 与 `_poly_add` 原始堆执行导出的精确乘积 `- 1` 长度界；
- C++ 使用的两个低半部长度 `min lenA m`、`min lenB m`。

`HgcdRecursiveLengthInvariant` 新增的正长度字段在真实 base/iterator
执行分支中分别由输入正长度和欧几里得循环不变量导出。首次重构定理显式消费
源入口的 `lenB < lenA`，没有把它藏进 workspace、规格 oracle 或额外公理。

## 验证

- `lake build CLPoly.Impl.StrictHGCDRawRefinement`
- `python3 proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `#print axioms`：严格下降算术定理仅依赖 `propext`、`Quot.sound`；
  原始堆精化定理仅再依赖 `Classical.choice`，无 `sorryAx`。

## 后续 leading-A 接口

原始表示长度桥同时推广为：只要整个低部的规范化长度严格小于
`shift + highLength`，移位高部的最高项就不可能被低部抵消，输出长度精确等于
`shift + highLength`。旧的 `lowLength ≤ shift` 版本保留为推论。这个推广
直接对应最终重构所需的次数分离，不引入“结果非零”之类的规格前提。

## 最终 A 首项保持

对照 C++ 第二次重构公式，新增并从真实迭代执行证明统一矩阵系数界：
四个 `_mat_one`/`_mat_row_update` 描述符长度均不超过
`inputLength - inputLength / 2`。每一步证明只使用源循环守卫、真实 divrem
商长度和两次真实行更新的长度结果。

递归长度契约现在同时携带 `aboveHalf` 和该系数界。于是最终重构的
`S[3] * lowA`、`S[1] * lowB` 均严格低于 `X^k * secondA` 的最高项，
由规范化 raw 表示推出 `result.lenA = k + second.lenA`，特别得到
`0 < result.lenA`。这不是对输出非零的假设，而是源矩阵更新与物理重构的推论。
代入源码的 `k = 2*m-lenB2+1` 后还得到
`outerLength / 2 < result.lenA`，因此该不变量可继续传给递归父调用。

## 组合矩阵的完整系数界

长度不变量进一步覆盖第 1、3 行与返回 A 的配对界；该界由真实
`_mat_row_update` 步骤逐次保持，而不是在递归出口假设。结合首次重构的精确
`lenA = m + first.lenA`，可以从真实 divrem 商界推出首次重构后的
`lenB ≤ lenA`。

最终的 quotient-update 先交换矩阵行，再分别更新两列。证明现在按列复用同一
算术引理，并同时闭合输出 0、1、2、3 四个描述符的半长界。其前提精确对应
`hgcdMatApplyQuotientEntries_length_bounds` 和
`hgcdRecursiveCombineMatrix_length_bounds` 给出的真实 C++ 描述符关系；没有用
规格矩阵、oracle、fuel 或 L2 结果替代生成代码的执行。
`hgcdRecursiveCombineMatrix_coeff_bounds` 已把四项算术界直接实例化到该真实执行
返回的 `modified` 与最终矩阵描述符，因而递归契约不再需要额外的规格侧矩阵界。

第一次成对重构现在还有执行级次序桥：同一次成功的四调用重构同时导出
`result.lenA = lenA / 2 + first.lenA`、A 正长度以及
`result.lenB ≤ result.lenA`。第二次递归调用因此可以直接消费真实描述符次序，
无需规格侧补充“输入已经排好序”的假设。

良基递归体的 proof-only provider 已相应加强为
`HgcdFirstReconstructionInvariant`。它一次携带精确 A 长度、A 正长度、真实操作数
次序和外层严格下降；`hgcdRecursiveBodyBelow` 的可执行分支仍完全相同，只从该
擦除结构读取下降证明。后续中间 divrem 精化可以直接读取同一个执行事实，而不
再另设 oracle 式前提。

## 最终尾部的操作数契约

真实 `hgcdRecursiveFinish` 成功执行现在直接导出最终 A/B 均不超过外层输入、
A 为正且严格高于外层半长、B 满足外层停止界。证明展开同一次
`hgcdRecursiveFinish_exec` 返回的第二重构，消费第二子调用的完整长度不变量和
源码 `k/c0` 公式；没有把最终长度作为 finish workspace 的假设。矩阵描述符部分
继续由实际 combine-matrix 执行定理提供，两部分将在总递归归纳步中合并。

## 完整 finish 长度不变量

最终矩阵界现已从“统一半长上界”加强为四个 sharp `row + finalA` 配对界。
证明保留真实商长度、第一子调用的 `row + firstA`、第二子调用的
`row + secondA`，以及源码恒等式 `k + c0 = reconstructedLenB`，并直接实例化
到 quotient-update 与矩阵乘法返回的描述符。

`hgcdRecursiveFinish_lengthInvariant` 随后在同一次 `hgcdRecursiveFinish_exec`
轨迹上合并操作数契约、四行 sharp 界和四项系数界，得到完整
`HgcdRecursiveLengthInvariant`。其中 `row1B/row3B` 由真实 B 停止界和 A
above-half 推出操作数次序后导出；没有额外假设最终 B 已小于 A。

## 良基 dispatch 归纳桥

`hgcdRecursiveDispatchBelow_rawInvariant` 已覆盖源码的 cutoff 分派：小输入臂展开
真实 iterator、stabilize 和 store 执行；大输入臂只消费严格更小 `lenA` 调用的
语义归纳假设。成功结果通过 proof-erased value 相等运输，不替换堆、矩阵或长度
字段。该桥把两个递归调用点统一为同一个无 fuel 归纳接口。

完整 `hgcdRecursiveBodyBelow` 的 base 路径也已闭合：从 body 的成功执行反演
真实 base helper，分别覆盖 `computeM=true/false` 的矩阵初始化差异，并将两次
原始复制的表示结果运输到 body 返回记录。该定理证明的是完整 body 分支，而不
只是孤立 base helper。

第一次重构命中 early guard 时，现在可直接构造完整父级长度不变量：最终 A 的
精确移位长度把第一子矩阵的四行配对界提升到外层输入，B/A 次序导出两条 row/B
界，而第一子调用的系数界单调提升为外层系数界。该结构来自真实重构不变量和
early guard，不作为 early workspace 的附加结果假设。

## 完整递归体的 early-return 路径

`hgcdRecursiveBodyBelow_early_rawInvariant` 现已把非 base 的 early-return
控制流闭合为共同 raw 不变量。第一次子调用的语义来自实际
`hgcdRecursiveDispatchBelow` 成功执行的良基归纳结果；随后用真实四调用成对重构
推出完整输入的 transform、带符号行列式和 GCD 保持，再用真实输出 memcpy 与可选
矩阵复制构造递归结果。

证明最终通过 `HgcdRecursiveResult.ext_value` 将 early helper 的返回记录运输到
完整 body 的实际返回值。物理 workspace 只描述有效性、分离性和前缀保持，不携带
预先选定的 L2 结果；执行路径中也没有计数参数或规格回退。

## 非 early 尾部的统一语义出口

新增 `hgcdRecursiveRawInvariant_of_finish_semantics`，把第一次递归矩阵、真实
middle divrem 等式、第二次递归矩阵和 finish 返回记录装配成共同 raw 不变量。
最终矩阵条目严格采用生成代码执行的
`first * ([[quotient,1],[1,0]] * second)` 顺序；transform、返回符号对应的带符号
行列式以及 GCD 保持均由该实际矩阵乘积推导。

该定理不创建输出或终止事实：A/B raw 表示、可选矩阵表示、停止界和完整长度
不变量都必须由实际 finish 执行及已有长度定理提供。它为下一步将
middle、第二次良基 dispatch 和 finish 控制流闭合到完整 body 提供统一出口。

## 定长低前缀的真实乘法语义

检查第二次重构输入时发现，源码在 `k` 处取得的低前缀允许末尾系数为零，因此不能
普遍满足“描述符长度已经 normalise”的 `RawDensePolyRep`。继续把该性质放进物理
workspace 会形成无法从真实执行导出的残留假设。

新增 `RawCanonicalPolySlice` 与 `hgcdRecursiveMulTerm_refines_slice`，直接调用乘法层
已有的 `mul_refines_slice`：输入只要求有效、machine-canonical 且完整可观察，输出
保留 C++ 的固定乘积容量和正确 L2 乘积，不提前声称结果已规范化。零长度分支也从
真实空切片推出零多项式。下一步将重构的乘积/加减链改用该接口，并只在实际
normalise 调用之后恢复 `RawDensePolyRep`。

## Middle 的完整低高分解

`hgcdRecursiveMiddle_split_reps` 现从同一次真实 divrem 返回堆导出第二递归与最终
重构需要的四个多项式表示。divisor/remainder 的低部使用允许尾零的
`RawCanonicalPolySlice`，`c0/d0` 高部则由原规范化多项式的 suffix 性质得到
`RawDensePolyRep`；两者分别满足精确的 `low + X^k * high` 分解。

当 `k` 超过 remainder 长度时，证明按源码得到零长度 `d0`，从物理有效空切片
推出零多项式，而不是补充一个规格侧零值。该桥消除了第二递归输入和 finish 低部
之间原先缺失的语义联系。

## Reconstruction 只在真实减法后规范化

`RawCanonicalPolySlice` 已下沉到公共 raw 多项式表示层，并新增从
`RawDensePolyRep` 的单向转换。`polySub_equalLength_refines`、左右长分支以及统一
`polySub_refines` 现在只要求两个输入是 canonical 定长 slice；其返回值仍由生成的
`dense_upoly_zp__poly_sub_ir` 内部真实 `normaliseU64` 得到
`RawDensePolyRep`。

`hgcdRecursiveReconstructB_refines` 和 `hgcdRecursiveReconstructA_refines` 已实际迁移
到 slice 乘法与 slice 减法接口。成对重构也允许低部为定长 slice，并在两个真实
乘积、符号选择减法和后续 lift-high 执行之后返回规范化 A/B。因此 final
reconstruction 不再依赖低前缀“碰巧已经规范化”的不可导出前提。

## 第二次重构的真实语义接口

最终重构的 operand、length 和 semantic 三层定理现统一接受
`RawCanonicalPolySlice` 低部；因此 `hgcdRecursiveMiddle_split_reps` 从真实 divrem
得到的 `b2[0:k]` 与 `d[0:k]` 可以直接进入 finish，不再要求规格侧补出
normalise 事实。`hgcdRecursiveReconstructPair_preserves_input` 与
`hgcdRecursiveFinish_refines` 的结论也加强为生成四调用实际构造出的精确多项式，
避免用无约束存在量掩盖第二次重构的操作数身份。

新增 `hgcdRecursiveFinish_preserves_input`：它反演同一次真实 finish 执行，使用
第二子调用的 raw 不变量和 middle 的低高分解，将高 suffix 上的 transform 提升为
完整 divisor/remainder 上的 transform，同时保留实际 A/B 输出、返回符号和可选
矩阵乘积。`hgcdRecursiveRawInvariant_of_finish_execution` 随后把这些具体字段装配为
共同递归不变量，为展开完整非 early body 只留下 middle/第二 dispatch 的物理
frame 接线。
