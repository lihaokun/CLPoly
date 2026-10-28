# SquarefreeZp 编译修复与 derivative_toPoly_eq 证明完成

## 日期
2026-05-30

## 做了什么
1. **完成 `derivative_toPoly_eq` 证明**：`SparsePolyZp.toPoly p (SparsePolyZp.derivative f) = Polynomial.derivative (SparsePolyZp.toPoly p f)` — 这是连接 C++ 和 L2 求导操作的核心引理。

2. **修复所有编译错误**：`Refinement/SquarefreeZp.lean` 从 ~30 个报错降至零报错（仅剩警告和 `sorry`）。

## 关键决策

### 1. 提取 `listSum_derivative_eq` 独立引理
将 List 层面的求导对应从 Array 层面分离，简化了证明。避免了 Array 归纳的复杂性和 `let` 变量 shadow 问题。

### 2. 使用 `List.filterMap_congr` 统一 prime
`SparsePolyZp.derivative` 使用 `f[0]!.snd.prime` 作为全局模数，而 `listSum_derivative_eq` 使用 `c.prime`（逐项）。通过 `UInt64.toNat_inj` 证明两者相等（因 `WellFormed` 保证所有 prime 的 `toNat = p`）。

### 3. 使用 `rw` 链避免 `calc` 的 `Zp.toZMod` 规范化问题
发现 `listSum_cons` 使用 `Zp.toZMod p c` 而非 `(c.val.toNat : ZMod p)`，导致 `calc` 因语法差异生成无法关闭的子目标。改用 `rw` 链解决。

### 4. 深度修复 UInt64→Nat→ZMod 的不一致
`derivative` 使用 UInt64 算术 (`c.val * m.deg.toUInt64 % p`)，而 `Polynomial.derivative` 在 ZMod p 中运算。需要 `h_no_overflow` + `h_deg_bound` 条件保证 UInt64 不溢出且 degree 不截断。

### 5. 主定理框架搭建
添加 `h_deg_bound` 参数，补完 `natDegree = 0` 分支，为 p-th root 分支和 Yun 分支预留结构。

### 6. 修复若干引理
- `normalization_toPoly_eq`: 使用 `listSum_filter_zero_coeff` 替代不存在的 `listSum_filter_zero`
- `listSum_filter_zero_coeff`: 处理 `Bool` vs `Prop` 的 `!=` / `≠` 转换
- `pair_vec_div` 调用补上 `()` 参数
- `upoly_make_monic` 调用补上 `.2` 提取
- `toPolyList_append`: 用 `++` 替代 `+`
- `toPolyList_scaleExp`: 标记为 `sorry`（原陈述数学上不正确，但未使用）

## 修改的文件
- `CLPoly/Refinement/SquarefreeZp.lean`：全部修改

## 剩余 `sorry`
- `get_deg_toPoly`（简单，未使用）
- `gcd_toPoly_eq`、`pair_vec_div_toPoly_eq`、`extract_pth_root_toPoly_eq`、`upoly_make_monic_snd_toPoly_eq`、`yunLoop_toPolyList_eq`（辅助引理）
- `toPolyList_scaleExp`（已废弃）
- `__squarefree_Zp_ir_refines` 主定理体（p-th root 分支和 Yun 分支）
