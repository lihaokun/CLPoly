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

## 下一步

继续闭合 `_mat_row_update` 的非零分支：接入已经完成的严格 `_mul` 和
`_poly_add` 精化，证明新矩阵行精确为 `old[i1] + Q*old[i0]`。之后才能证明
`_hgcd_iter` 每轮同时保持矩阵变换关系和 GCD。
