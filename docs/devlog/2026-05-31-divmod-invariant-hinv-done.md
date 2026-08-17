# divmodAux 代数不变量证明（h_inv_next）完成

## 日期
2026-05-31

## 做了什么
1. **新增 `mem_getFirst_toList` 引理**：非空 SparsePolyZp 的首元素 `r[0]!` 在 `r.toList` 中。使用 `getElem!_def` + `getElem?_def` + `hpos` 处理 `r[0]!` vs `r[0]` 的差异，避免了 `Array.get` 不可用的问题。
2. **完成 `h_coeff_red` 证明**：`coeff = r[0]!.snd * lc_g_inv` 是 `Zp.Reduced p`：
   - `coeff.prime.toNat = p`：通过 Zp 乘法的 prime 来自左操作数的定义事实（`rfl` 可证）
   - `coeff.val.toNat < p`：通过 `Nat.mod_lt` + `UInt64.toNat` round-trip 性质
3. **完成 `h_inv_next` 代数不变量**：使用 `hq'`（`listSum_append` 证明 `toPoly q' = toPoly q + toPoly term`）、`hsub`（`toPoly_sub`）、`hmul`（`toPoly_mul`），通过 `calc` 证明 `f = q'*g + r'`。
4. **填充 `hred_q'` 和 `hred_r'`**：
   - `hred_q'`：`q.push term` 保持 AllReduced（`q.toList ++ [term]` 中每个元素都 Reduced）
   - `hred_r'`：由 `WellFormed_arr.sub` 直接给出
5. **重构 `h_aux` 签名**：增加 `hred_q` 和 `hred_r` 参数，支持 toPoly homomorphism lemmas。
6. **修改 `divmodAux_invariant` 签名**：增加 `hred_q` 和 `hred_r` 参数。
7. **新增 `allReduced_empty` 引理**：空数组总是 AllReduced。

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 当前状态
- 构建零错误
- `divmodAux_invariant` 中的三个 `admit` 已全部填充：
  - ✅ `h_inv_next`：填充
  - ✅ `hred_q'` 和 `hred_r'`：填充
  - ❌ `h_deg_dec`：留为 `admit`（需 leading term cancellation 证明）
- `divmodAux'_wf` `decreasing_by`、`divmodAux'_wf_eq`、`gcdAux'_wf` `decreasing_by` 仍为 `admit`
- `divmod_deg_lt`、`divmod_snd_wellFormed`、`divmod_snd_allReduced`、`gcdAux_unfold` 仍为 `admit`
