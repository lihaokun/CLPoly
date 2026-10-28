# SquarefreeZp 主定理骨架 + 关键发现

## 日期
2026-05-30

## 做了什么
1. `extract_pth_root_toPoly_eq` 框架证明 — 使用 `expand_contract` + `expand_injective` 的往返法
2. 补充 `gcd_toPoly_eq` 签名（含 `h_p2` 条件）
3. 发现 C++ `__squarefree_Zp_ir_def` 对常数多项式不终止

## 关键发现

### C++ 模型对常数多项式不终止
`__squarefree_Zp_ir_def` 对常数多项式的递归不终止，因为 p-th root 分支 (`derivative = 0`) 下：
- `extract_pth_root` 把 `deg=0` 映射为 `deg=0`（`0 / p = 0`）
- `upoly_make_monic` 使首项为 1（对常数多项式就是它自己）
- 递归调用 `__squarefree_Zp_ir_def` 参数不变 → 死循环

`decreasing_by` 中的两个 `sorry` 确认了终止性未证明。

这意味着 `natDegree = 0` 分支的 refinement 在现有模型下**不可证明**。修复方案：
- 在 `__squarefree_Zp_ir_def` 开头加 `if get_deg f = 0 then #[]` 检查
- 或证明主定理的输入永远不是常数（通过调用方约束）

### extract_pth_root 的证明策略
用 `expand p` 做往返：
```
expand p (toPoly (extract f)) = toPoly f    -- 待证
expand p (contract p (toPoly f)) = toPoly f  -- 由 expand_contract
expand 单射 → toPoly (extract f) = contract p (toPoly f)
```
`h_expand_L` 需要 extract 循环的列表变换引理。

## 当前状态
- ✅ 全文件零编译报错
- ❌ 6 个 `sorry`（`gcd_toPoly_eq`、`extract_pth_root_toPoly_eq`、`yunLoop_toPolyList_eq`、主定理 3 分支）
