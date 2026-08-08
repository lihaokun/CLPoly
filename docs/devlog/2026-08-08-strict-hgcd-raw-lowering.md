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
- 四项 `HgcdMatPolyRep` 现可由逐项 layout/prefix frame 整体传递，并已
  专门连接到真实 recursive memcpy；矩阵切片有效性保持为显式 L1 条件。

下一步证明初始化 copy 的表示传递，以及每轮 divrem/两次行更新共同保持
Euclid 对与矩阵变换关系，最终导出循环停止条件和 GCD 不变式。
