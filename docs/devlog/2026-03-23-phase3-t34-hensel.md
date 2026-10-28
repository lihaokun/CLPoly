# Phase 3 T3.4 Hensel 提升完成 — 0 sorry

> 日期：2026-03-23
> 分支：`feature/formal-proofs`

---

## 做了什么

完成 Hensel 提升的 Lean 4 形式化证明（366 行，0 sorry）。

**证明内容**：
- `hensel_step`：2-factor 单步 mod m → mod m²（ℤ[x] 方法）
- `zmod_ker_mul_eq_zero`：ZMod(m²) 中 ker(π)² = 0
- `poly_ker_mul_eq_zero`：多项式级 ker² = 0
- `isCoprime_lift_sq`：IsCoprime 从 mod m 传播到 mod m²
- `hensel_two_factor`：迭代 mod p → mod p^k（对 k 归纳 + 投射）

**核心数学**：
- ℤ[x] 方法：g' = g + C(m)·τ, h' = h + C(m)·σ
- 恒等式 `g'·h' - f = C(m²)·(something)` 用 `ring` + `linear_combination` 验证
- C(m²) 在 ZMod(m²) 中映射为 0 → 乘积保持
- IsCoprime 传播：Bézout lift 给 1+δ，δ²=0 → nilpotent → unit

## 为什么做

ZZ 因式分解管线的关键组件。HenselCorrect spec 要求将 mod p 因式分解提升到 mod p^k。

## 关键决策

1. **ℤ[x] 方法而非 ZMod(m²) 内部方法**：在 ℤ[x] 中构造 g', h'，避免在 ZMod(m²) 中提取 `∃ a, x = m*a`（Lean 中很困难）。关键恒等式 `g'*h' - f = C(m²)*(...)` 完全在 ℤ[x] 中用 `ring` 验证。

2. **IsCoprime 传播两层方案**：
   - mod m → mod m²：ker(π)²=0 → δ²=0 → `IsNilpotent.isUnit_one_add`
   - mod p → mod p^k：每个系数是 p 的倍数 → nilpotent（`p^k=0`）→ `Polynomial.isNilpotent_iff` → `IsNilpotent.isUnit_one_add`

3. **迭代用 k 归纳（非 doubling）**：hensel_step(m=p^k) 给出 mod p^{2k}，投射到 p^{k+1}（因 k+1 ≤ 2k when k ≥ 1）。

4. **不含度数保持**：当前 hensel_two_factor 只证乘积还原 + mod p 保持。度数保持（HenselCorrect 条件 3）需要 divByMonic 控制，留作后续。

## 遇到的问题

| 问题 | 解法 |
|------|------|
| `linarith` 对多项式无效 | `sub_eq_iff_eq_add'.mp` 或 `eq_add_of_sub_eq` |
| `ring` 无法用 hypothesis | `linear_combination C m * e_int * hbez_eq` |
| `ZMod.natCast_zmod_eq_zero_iff_dvd` 不存在 | 正确名 `ZMod.natCast_eq_zero_iff`（ℕ版）、`intCast_zmod_eq_zero_iff_dvd`（ℤ版） |
| `Int.ediv_mul_cancel` 方向 | `Int.mul_ediv_cancel'` 给 `a * (b/a) = b`（+ `.symm`） |
| `ZMod.cast_eq_val` + `natCast_zmod_val` 转换 | 显式 `rw` 链 |
| `Polynomial.map_surjective` 需显式 ring hom | 提供 `(Int.castRingHom (ZMod m))` |
| `Nat.one_lt_pow` 参数 | 接受 `k ≠ 0` 不是 `0 < k` |
| `p^(2*k)` vs `(p^k)^2` | `hdvd_pk2 : p^(2*k) = (p^k)^2` by `ring` |
| `Polynomial.isNilpotent_iff` | 逐系数证 nilpotent |
| nl-proof 4 轮审核 | ker²=0 + divByMonic 度数 + ℤ[x] 翻译方案 |

## 涉及的文件

| 文件 | 状态 |
|------|------|
| `proof/lean/CLPoly/Algorithm/Hensel.lean` | **新建**，366 行，0 sorry |
| `proof/lean/CLPoly.lean` | 添加 import |
| `proof/nl-proof/phase3-t34-hensel.md` | nl-proof v4（4 轮审核） |
| `proof/CLAUDE.md` | 添加 nl-proof 审核标准 |

## 度量

- 耗时：~8 小时（nl-proof 4 轮审核 ~4h + 形式化 ~4h）
- 迭代：nl-proof 4 版本 + Lean ~15 轮编译修复
- Lean 新增行数：366 行
- 放弃的方案：
  - ZMod(m²) 内部方法（`∃ a, x = m*a` 提取困难）→ 改用 ℤ[x] 方法
  - 后验度数论证（`natDegree_map_of_leadingCoeff_ne_zero` 方向错误）→ 暂不含度数保持

## 最终状态

| 模块 | 行数 | sorry |
|------|------|-------|
| SQF | 1501 | 0 |
| DDF | 435 | 0 |
| EDF | 257 | 0 |
| Hensel | 366 | 0 |
| Pipeline | 379 | 0 |
| **合计** | **2938** | **0** |

`lake build` 1913 jobs 全通过。
