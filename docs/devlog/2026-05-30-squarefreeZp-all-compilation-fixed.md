# SquarefreeZp.lean 全文件编译通过 + 主定理结构搭建

## 日期
2026-05-30

## 做了什么
1. 修复 `Refinement/SquarefreeZp.lean` 全部编译错误（从 ~30 个降至零）
2. 修正了多个引理表述（`upoly_make_monic_snd_toPoly_eq`、`pair_vec_div_toPoly_eq`）
3. 利用 `ZpArith.lean` 中已有定理（`__upoly_make_monic_ir_refines`、`__upoly_divmod_ir_refines`）填充引理
4. 搭建了主定理 `__squarefree_Zp_ir_refines` 的证明骨架

## 关键修复

### 表述修正
- `pair_vec_div_toPoly_eq`：原说 `= toPoly p q`，实际 `pair_vec_div` 忽略 `q`，结果是 `(divmod f g).fst`
- `upoly_make_monic_snd_toPoly_eq`：原说 `= toPoly p f`，实际是 `= toPoly p (makeMonic f)`（首一化，非恒等）

### 利用已有定理
- `__upoly_make_monic_ir_refines` → 填充 `upoly_make_monic_snd_toPoly_eq`
- `__upoly_divmod_ir_refines` → 填充 `pair_vec_div_toPoly_eq`

### 识别 UInt64 != vs Prop ≠ 问题
`simp` 不能直接用 `hc0 : c.val ≠ 0`（Prop）重写 `c.val != 0`（Bool），需要 `simpa` 或额外转换。

### 识别 `calc` 的 `Zp.toZMod` 规范化问题
`listSum_cons` 使用 `Zp.toZMod p c` 而非 `(c.val.toNat : ZMod p)`，导致 `calc` 生成不可关闭子目标。改用 `rw` 链解决。

## 剩余 `sorry`（7个）
| # | 引理 | 难度 | 依赖 |
|---|------|------|------|
| 1 | `get_deg_toPoly` | 中 | 不用于主定理 |
| 2 | `gcd_toPoly_eq` | 高 | 需证明 `toPoly` 保持欧几里得算法 |
| 3 | `extract_pth_root_toPoly_eq` | 高 | extract 循环语义对应 contract p |
| 4 | `yunLoop_toPolyList_eq` | 很高 | Yun 循环对应 |
| 5 | `toPolyList_scaleExp` | 表述有误 | 已废弃，不用于主定理 |
| 6 | 主定理 `natDegree=0` 分支 | 中 | C++ 对常数多项式返回 `#[]` |
| 7 | 主定理 `derivative=0` 和 `derivative≠0` 分支 | 高 | 依赖引理 2-4 |

## 下次方向建议
1. 最难但最关键：`gcd_toPoly_eq` — 需要 `toPoly` 的欧几里得算法保持性
2. 或 `extract_pth_root_toPoly_eq` — extract 循环语义
3. 或直接填主定理的 `natDegree=0` 分支
