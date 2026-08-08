# HGCD 变换不变式单步

日期：2026-08-08

## 做了什么

- 定义 C++ HGCD 矩阵的精确方向：原始多项式对等于矩阵乘以当前 Euclid 多项式对。
- 证明一次除法恒等式和源码两次行更新保持该变换关系。
- 证明两次源码行更新使 2×2 行列式符号翻转。
- 加强零分支执行定理，给出交换后两个矩阵描述符的精确指针和长度。
- 闭合零分支语义：实际描述符交换表示 `[entry1 + quotient*entry0, entry0]`。
- 统一零分支与非零 `_mul`/`_poly_add` 分支，得到覆盖生成函数全部控制流的 `matRowUpdate_refines`。
- 将 identity 四个矩阵条目加强为规范化 `RawDensePolyRep`，并证明两个实际初始化 `memcpy` 保持该四条目不变式。
- 将 identity raw 表示与代数方向合并，证明初始化状态满足“原始对 = identity × 当前对”。
- 证明完整 row update 的通用内存帧：与 `T`、`scratch`、`t` 写区分离的 raw 多项式在零/非零全部执行分支中保持不变。
- 将 C++ `sgn` 与矩阵行列式 `±1` 关联，证明行更新和符号取反同步保持，并由带符号变换推出规范化 gcd 保持。
- 将非终止循环执行分解与 `polyDivrem_next_state` 合并，绑定同一个实际 divrem 返回值的 quotient/remainder raw 表示、gcd 保持、两次 row update 和递归尾调用。
- 证明 row update 未选中索引的指针/长度描述符保持，并组合两次更新得到 `(0,1)` 与 `(2,3)` 的交叉描述符帧。
- 定义四条目的 `hgcdStepEntries`，证明其同时保持 HGCD 变换关系和同步更新 signed determinant。
- 将每次 row update 所需的缓冲区有效性和别名条件封装为纯 L1 的 `MatRowUpdateWorkspace`，不在其中预设任何 L2 运算结果。
- 证明两个实际生成的 `_mat_row_update` 调用共同得到完整四条目 `hgcdStepEntries` raw 表示；证明过程中显式保持 quotient、未选中矩阵条目及描述符，而不是直接改写为数学矩阵公式。
- 加强共享的 `polyDivrem_next_state` raw→safe 桥，保留真实 `_poly_divrem` 已证明的除法恒等式，并将它传入 HGCD 单轮接口；此前该接口只保留 gcd 等式，无法闭合矩阵变换不变量。
- 定义外部 live polynomial 的 `MatRowUpdateGuardWorkspace`，并证明任意满足物理分离条件的 raw 多项式经过两次实际 row update 后保持；这将用于把当前 divisor 和 divrem remainder 传给下一轮 HGCD 状态。
- 将 divrem 的递减 quotient recursion 内存帧从专用 divisor 推广到任意与 `Q/W3` allocation 分离的 UInt64 slice；证明逐次跟随实际 `Q` 写入、可选 `addMulLoop` 的 `W3` 写入及递归调用，为矩阵跨 divrem 保持奠定基础。
- 闭合完整 `_poly_divrem` 的任意外部 slice 帧，覆盖短输入 `copyU64` 与长输入三段循环，并据此证明四项 HGCD raw 矩阵跨真实 divrem 保持。
- 将实际 divrem、两次实际 row update、raw 状态轮换、矩阵变换、signed determinant、规范化 gcd 及严格 `lenR < lenB` 组合成 `hgcdIterationCalls_refine` 单轮精化定理。
- 定义不含 L2 结果的循环物理工作区 provider，并以生成函数相同的 `state.lenB` 度量证明完整 `hgcdIterLoop_refines`：每轮调用上述真实单轮定理，以实际 `lenR < lenB` 进入递归，最终得到停止条件、最终 raw 状态、transform、signed determinant 和全程 gcd 保持。
- 将真实 `_mat_one`、两次有序 `copyU64` 初始化与良基循环组合为完整 `hgcdIter_refines`，从输入 raw 多项式直接得到生成 `_hgcd_iter` 最终状态的不变量、gcd 保持和源码停止界。
- 开始 `_hgcd_recursive` 严格 lowering：提取并实现真实 base branch 的可选 `_mat_one` 与按源码顺序执行的 `A ← a`、`B ← b` raw copies；证明 compute-matrix 分支复用同一初始化执行并得到 identity 矩阵及两个输出 raw 表示。
- 闭合 base branch 的 `compute_M=false` 路径：不调用任何矩阵规格或实现，直接跟随两次 raw copy，并用第一次 copy 对 `b` 的帧和第二次 copy 对 `A` 的帧证明源码声明的 `{B,a}` 别名安全顺序。
- 严格降低非 base 分支的工作区切片：记录 `a2/b2/a3/b3/q/d/T0/T1`、R/S 八个矩阵条目、三肢 `W3` 和 `W_next` 的源码指针推进，并用 `hgcdRecursiveWorkspace_layout` 固定每个等式及两个矩阵描述符的四项布局。
- 降低第一次递归/iterator 的高半部输入 `a0=a+m`、`b0=b+m` 及源码条件长度，并从非 base 条件与 `len_b < len_a` 证明 `len_a0 < len_a`，为 `_hgcd_recursive` 的真实良基递归提供下降证据。
- 将 iterator 返回后的 R/S 稳定化块降低为两个以 `4-i` 为度量的 raw 循环：第一遍逐项 `copyU64` 到连续 stage，第二遍逐项复制回调用前指针并同步恢复矩阵 descriptor；组合函数保持两个循环及其错误传播顺序。
- 证明第二个真实恢复循环的每个成功后缀都保持矩阵描述符有效性，并逐字保持 iterator 返回的四项长度；该结论由实际四次 copy/descriptor 更新的良基循环归纳得到，不预设 L2 矩阵结果。
- 提取恢复循环的纯 descriptor 效果并证明其保持数组长度；再证明每个成功 raw 执行返回的指针数组精确等于这一逐项 `Array.set` 效果，为最终四指针复原与 staged raw 语义搬运提供执行桥。
- 展开起始索引 `0` 的四次实际 descriptor 写入，并逐索引证明最终指针数组等于调用前保存数组；组合两个真实 copy 循环，得到完整 stabilization 成功时“旧指针全部恢复、iterator 新长度全部保留、descriptor 仍有效”的执行定理。
- 固定 stage 缓冲区中四项的源码顺序偏移及总长度，并证明真实第一循环的最终 offset 恰为四项长度之和。
- 定义只含容量和别名条件的 `HgcdMatStabilizeWorkspace`；以生成循环的 `4-i` 良基度量逐次调用真实 `copyU64`，证明第一循环成功终止、保留当前矩阵 raw 表示，并在连续 stage 切片上建立全部四个规范化多项式表示。既有 stage 项通过相邻切片不重叠逐轮保持，workspace 不携带任何 L2 结果。
- 对第二个生成 restore 循环建立同样以 `4-i` 递减的 raw 语义归纳：每轮从对应 stage 切片复制到保存的原始指针，保持整个 stage 和此前恢复的目标项，并同步跟随真实 descriptor 更新。
- 从纯物理 workspace 单独证明 restore 循环必然成功终止，再与 descriptor 恢复定理组合，得到完整 `hgcdMatStabilize` 的总精化：两段实际循环执行成功，输出矩阵在恢复后的旧指针和 iterator 产生的新长度上表示同一组四个 L2 多项式。
- 继续降低 iterator 后的交叉保护块：当 `pA != a3 && pB == a3` 时严格先执行 `b3 ← pB` 再执行 `a3 ← pA`，其余情况按源码分别跳过或执行两次条件 copy；两个执行形状定理固定了真实分支和错误传播顺序。
- 闭合交叉保护块的 raw 多项式语义：危险分支在第一次写 `b3` 时显式 framing `pA`，第二次写 `a3` 时再 framing 已保存的 `b3`；普通分支分别覆盖两次 copy、跳过 A、跳过 B、两次均跳过，并统一推出最终 `a3/b3` 表示 iterator 返回的两个多项式。

## 为什么做

HGCD 循环的 raw 执行只有在每轮都保持“原始对—矩阵—当前对”的关系时，才能最终连接到 GCD 及 SQF。零分支也必须从真实 L1 表示推导，不能当作特殊规格捷径跳过。

## 关键决策

- 采用源码实际的右乘方向：每一行 `[u,v]` 更新为 `[v+Q*u,u]`，吸收 `A=Q*B+R`。
- 零分支中的消失乘积由零长度 `SlicePolyRep` 推出 quotient 或 entry 为零，而不是增加 L2 前提。
- 行列式翻转单独证明，后续将与状态中的 `sgn := -sgn` 对齐。
- 两次行更新之间复用同一 quotient；因此把 quotient 的内存帧作为独立定理证明，并要求 addition destination 与 quotient allocation 分离。
- 循环级接口明确分离语义不变量与物理工作区：provider 只量化成功 L1 调用并给出有效性/容量/别名事实，不携带 quotient、remainder、矩阵公式或其他 L2 结论。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictHGCDRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-08-hgcd-transform-step.md`
