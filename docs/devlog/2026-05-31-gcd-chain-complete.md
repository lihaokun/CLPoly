# extGcdAux 证明链完成 + admits 清理

## 日期
2026-05-31

## 做了什么
1. **完成 `extGcdAux_linearity_gcd`**：扩展了线性性引理，增加了 `g ∣ old_r` 和 `g ∣ r` 的整除性质。
   证明用强归纳于 `r.natAbs`，是目前最完整的 extGcdAux 性质引理。

2. **完成 `extGcdAux_bezout`**：由 `extGcdAux_linearity_gcd` 直接推出 `coeff * B + u * A = g`。

3. **添加了骨架引理**：
   - `extGcdAux_gcd_nonneg`：admit（非负输入返回非负 gcd）
   - `Zp_toZMod_inv_mul_self`：admit（Zp 逆的正确性）
   - `divmod_deg_decrease`：admit（多项式长除法的度递减性质）

4. **清理了大量遗留代码**：删除了旧的 `extGcdAux_linear_full` 引理和不完整证明尝试。

## 当前证明链
```
extGcdAux_linearity_gcd (done)
    → extGcdAux_bezout (done)
    → extGcdAux_gcd_nonneg (admit) + Zp_toZMod_inv_mul_self (admit)
        → divmod_deg_decrease (admit)
            → h_deg_dec 在 divmodAux_invariant (done)
                → divmod_identity (需要 divmodAux'_wf_eq admit)
```

## 剩余关键 admit
| admit | 位置 | 阻塞 |
|-------|------|------|
| `extGcdAux_gcd_nonneg` | 近顶部 | extGcdAux 返回非负 gcd |
| `Zp_toZMod_inv_mul_self` | 近顶部 | Zp 逆正确性 |
| `divmod_deg_decrease` | 近顶部 | 长除法度递减 |
| `divmodAux'_wf_eq` | ~431 | _wf 与 partial def 等价 |
| `divmod_deg_lt` | ~666 | 余式度 < 除式度 |
| `divmod_snd_wellFormed` | ~672 | 余式保持 WellFormed |
| `divmod_snd_allReduced` | ~680 | 余式保持 AllReduced |
| `gcdAux_unfold` | ~686 | gcdAux 展开引理 |
| 其他 | ~690+ | gcd 性质、extract_pth_root 等 |

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 当前状态
- 构建零错误
- `extGcdAux_linearity_gcd` 和 `extGcdAux_bezout` 是本次完成的关键引理
- 剩余 admit 约 15 个
