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

## 下一步

继续闭合 `_mat_row_update` 的非零分支：接入已经完成的严格 `_mul` 和
`_poly_add` 精化，证明新矩阵行精确为 `old[i1] + Q*old[i0]`。之后才能证明
`_hgcd_iter` 每轮同时保持矩阵变换关系和 GCD。
