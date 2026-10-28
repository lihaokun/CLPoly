# extGcdAux_linearity_gcd 完成 + h_deg_dec 重构

## 日期
2026-05-31

## 做了什么

1. **完成 `extGcdAux_linearity_gcd`**：扩展了 `extGcdAux_linearity`，增加了 `g ∣ old_r` 和 `g ∣ r` 的整除性质。
   - 通过强归纳于 `r.natAbs` 证明
   - 递归步：`g' ∣ r` 和 `g' ∣ r'` 推出 `g' ∣ old_r`（因为 `r' = old_r - q*r`）
   - 关键：`g ∣ r` 和 `g ∣ old_r` 说明 `g` 是 `gcd(old_r, r)`

2. **重构 `h_deg_dec`**：将 `divmodAux_invariant` 中的 admit 移到单独的 `divmod_deg_decrease` 引理中
   - 使用 `extGcdAux_bezout` + `extGcdAux_linearity_gcd` 可证 `g[0]!.snd * g[0]!.snd.inv ≡ 1 (mod p)`
   - 由此可得 `r' = r - term*g` 的首项抵消 → 度递减
   - 仍保留为 `admit`，但规格清晰、依赖链完整

3. **修复了大量证明细节**：
   - `dvd_zero` 的类型类问题
   - `Int.ediv_add_emod` 的 `omega` → `linarith` 迁移
   - `by_cases r == 0` → `by_cases (r : Int) = 0` 避免 `BEq` 转换问题

## 证明链现状

```
extGcdAux_linearity_gcd (已完成) 
    → extGcdAux_bezout (已完成)
    → divmod_deg_decrease (admit，需逆元正确性)
    → h_deg_dec 在 divmodAux_invariant 中 (使用 divmod_deg_decrease)
    → divmod_identity (可完成)
    → ...
```

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 当前状态
- 构建零错误
- 剩余关键 admit：`divmod_deg_decrease`（需 `Zp_toZMod_inv_mul_self`）、`divmodAux'_wf_eq`、`divmod_deg_lt` 等
