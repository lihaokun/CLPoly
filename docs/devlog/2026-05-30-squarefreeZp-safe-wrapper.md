# SquarefreeZp safe wrapper + 状态整理

## 日期
2026-05-31

## 做了什么
1. 添加 `__squarefree_Zp_ir_safe` 安全包装解决 C++ 模型对常数多项式不终止的问题
2. 主定理改用 safe wrapper，`natDegree=0` 分支可直接用 `simp` 证明
3. 剩余 5 个 `sorry` 全部退回干净状态

## safe wrapper 方案
```lean4
noncomputable def __squarefree_Zp_ir_safe (p : ℕ) [Fact (Nat.Prime p)] (f : SparsePolyZp) : 
    Array (SparsePolyZp × UInt64) :=
  if (SparsePolyZp.toPoly p f).natDegree = 0 then #[]
  else Generated.__squarefree_Zp_ir f
```

- `natDegree = 0` 时直接返回 `#[]`（与 `sqfZp` 一致）
- `natDegree > 0` 时委托给 C++ 模型（此时不会进入常数死循环）
- 主定理 `__squarefree_Zp_ir_refines` 改为证明 safe wrapper 版本的 refinement

## 剩余 5 个 `sorry`

| # | 引理 | 难度 | 所需的关键子引理 |
|---|------|------|-----------------|
| 1 | `gcd_toPoly_eq` | 高 | `toPoly_mul` + 欧几里得算法保持性 |
| 2 | `extract_pth_root_toPoly_eq` | 高 | extract 循环的列表变换语义 + `(deg/p)*p = deg` |
| 3 | `yunLoop_toPolyList_eq` | 很高 | Yun 循环对应（当前表述可能有误） |
| 4 | 主定理 `derivative=0` 分支 | 中高 | 依赖 2 + `upoly_make_monic_snd_toPoly_eq`（已有） |
| 5 | 主定理 `derivative≠0` 分支 | 高 | 依赖 1 + 3 + `pair_vec_div_toPoly_eq`（已有） |

## 下次建议
`gcd_toPoly_eq` 是最基础的代数引理，建议先攻。需要 `h_p2` 条件和 `toPoly_mul`。
