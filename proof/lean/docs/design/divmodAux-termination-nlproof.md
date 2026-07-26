# divmodAux 终止性 nl-proof 草稿

日期：2026-07-26
目标：填掉 `CLPoly/Model.lean:1141` 的 `decreasing_by admit`，使 `divmodAux` 成为真 WF 递归（0 admit）。

## 0. 背景与度量

```lean
def divmodAux (g dg lc_g_inv) (h_lc : (lc_g_inv * g[0]!.snd).val = 1) (h_sorted_g)
    (q r) (h_sorted_r : Sorted r) : ... :=
  if hr : r.isEmpty then (q, r)
  else
    let dr := r[0]!.fst.deg
    if hd : dr < dg then (q, r)
    else
      let coeff := r[0]!.snd * lc_g_inv
      let d := dr - dg
      let term := #[(⟨d⟩, coeff)]
      let r' := r - (term * g)         -- = subImpl r (term*g)
      let q' := q.push (⟨d⟩, coeff)
      divmodAux ... q' r' ...
termination_by (if r.isEmpty then 0 else r[0]!.fst.deg + 1)
```

递归调用只在 `¬ r.isEmpty ∧ ¬ dr < dg`（即 `dg ≤ dr`）时到达。此时
`measure r = dr + 1`（`dr = r[0]!.fst.deg`）。要证

```
measure r' = (if r'.isEmpty then 0 else r'[0]!.fst.deg + 1) < dr + 1.
```

- `r'` 空：`0 < dr + 1` ✓。
- `r'` 非空 sorted：只需 `∀ z ∈ r'.toList, z.fst.deg < dr`（头是成员，得 `r'[0]!.fst.deg < dr`，即 `+1 < dr+1`）。

**核心命题 (ALL-DROP)**：`∀ z ∈ (subImpl r (term*g)).toList, z.fst.deg < dr`。

## 1. 首项抵消需要的不变量

`subImpl r (term*g)).toList = mergeAdd r.toList (negImpl (term*g)).toList`。
`r` sorted 非空 ⇒ `r.toList = a :: rest_r`，`a = r[0]!`，`a.fst.deg = dr`，`rest_r` 全 `< dr`。

要 mergeAdd 把首项消掉，需 `(term*g)` 的首项 `= (⟨dr⟩, c)` 且 `c.val = a.snd.val`，其余 `< dr`。
这需要 **两个新不变量**（三塔已备好保持引理）进 `divmodAux` 签名：

- `pm : UInt64`、`hq : 0 < pm.toNat`
- `ReducedB r pm`（每系数 `val < pm ∧ prime = pm`）、`ReducedB g pm`
- `NonZeroB r`（每系数 `val ≠ 0`）、`NonZeroB g`
- `h_lc_prime : lc_g_inv.prime = pm`
- `h_dg : g[0]!.fst.deg = dg`、`hg_ne : ¬ g.isEmpty`

递归调用处需证这些对 `r'` 保持：
- `ReducedB r' pm`：`ReducedB_subImpl` + `ReducedB_single_mul`（需 `coeff.prime = pm`）。
- `NonZeroB r'`：`NonZeroB_subImpl` + `NonZeroB_single_mul`（需 `ReducedB (term*g) pm`）。
- `Sorted r'`：已有（`subImpl_sorted + single_mul_sorted`）。
- `pm, g, dg` 相关的都不变（`g`、`dg`、`pm` 递归中恒定）。

`coeff.prime = pm`：`coeff = a.snd * lc_g_inv`，`Mul` 取左操作数 prime，故 `coeff.prime = a.snd.prime = pm`（`ReducedB r` 头）。

## 2. term*g 的首项（引理 A）

`term * g = addImpl #[] (scaleByMonomial ⟨d⟩ coeff g)`（`rfl`，`mulImpl` 单元素折叠）。
`.toList` 与 `(scaleByMonomial ⟨d⟩ coeff g).toList` 相同（`mergeAdd [] xs = xs` + toList/toArray 往返）。

`coeff.val ≠ 0`（下面证）⇒ `scaleByMonomial ⟨d⟩ coeff g = g.filterMap Gs`，
`Gs (mf,cf) = if (coeff*cf).val = 0 then none else some (⟨d + mf.deg⟩, coeff*cf)`。

`g` 非空 sorted ⇒ `g.toList = gh :: gt`，`gh = g[0]!`，`gh.fst.deg = dg`（`h_dg`），`gt` 全 `< dg`。

- 头 `Gs gh = some (⟨d + dg⟩, coeff*gh.snd)`。`d + dg = dr`（`d = dr-dg, dg ≤ dr`）。
  `(coeff*gh.snd).val = a.snd.val`（引理 B）`≠ 0`（`NonZeroB r` 头）⇒ 不被 filter 掉，确为头。
- 尾 `gt.filterMap Gs`：`sortedListB_filterMapShift` 第一部分，shift = d，输入界 dg，输出全 `< d + dg = dr`。

故 `(term*g).toList = (⟨dr⟩, c) :: tail`，`c.val = a.snd.val`，`c.prime = pm`（`Zp_mul_reduced`：`.prime` 取左 = coeff.prime = pm），`tail` 全 `< dr`。

## 3. Zp 模算术（引理 B、C）

**引理 B `Zp_mul3_lc`**：`a lc gh : Zp`，`a.prime = lc.prime = gh.prime = pm`，`a.val < pm`，
`(lc*gh).val = 1` ⇒ `((a*lc)*gh).val = a.val`。
证：设 `P = pm.toNat`。逐层 `Nat.mul_mod` + `mul_assoc` + `Nat.mod_mod`：
`((a*lc)*gh).val.toNat = (a.val.toNat * lc.val.toNat * gh.val.toNat) % P`
`= (a.val.toNat * ((lc.val.toNat*gh.val.toNat)%P)) % P`（mul_assoc + mul_mod）
`= (a.val.toNat * 1) % P`（`(lc*gh).val = 1 ⇒ (lc.val*gh.val)%P = 1`）
`= a.val.toNat % P = a.val.toNat`（`a.val < pm`）。toNat 单射 ⇒ `.val` 相等。

**引理 C `Zp_add_self_neg_val`**：`a c : Zp`，`c.val = a.val`，`c.prime = a.prime`，`a.val < a.prime`，
`0 < a.prime.toNat` ⇒ `(a + (-c)).val = 0`。
证：`(-c).val.toNat = ((a.prime - a.val) % a.prime).toNat`（`c=a` 代入）。`a.val < a.prime` ⇒ UInt64 减不回绕
`(a.prime - a.val).toNat = P - a.val.toNat`。`% P`：`a.val.toNat=0 → 0`；否则 `P - a.val.toNat`。
`(a.val.toNat + (-c).val.toNat) % P`：两情形都 `= P % P = 0`（omega 分情形）。`.toUInt64 0 = 0`。

## 4. mergeAdd 首项抵消（引理 D）

**引理 D `mergeAdd_cancel_lead`**：`a b`，`a.fst.deg = b.fst.deg = d`，`(a.snd+b.snd).val = 0`，
`rest_a rest_b` 全 `< d` ⇒ `∀ z ∈ mergeAdd (a::rest_a) (b::rest_b), z.fst.deg < d`。
证：`a.deg = b.deg` ⇒ mergeAdd 走 else-else 分支，`s = a.snd+b.snd`，`s.val=0` ⇒ `= mergeAdd rest_a rest_b`。
`mergeAdd_lt_all d rest_a rest_b` 给全 `< d`。

## 5. 组装 ALL-DROP

- `a := r[0]!`，`a.fst.deg = dr`，`rest_r` 全 `< dr`（Sorted r + sortedListB_iff）。
- 引理 A：`(term*g).toList = (⟨dr⟩, c) :: tail`，`c.val = a.snd.val`，`c.prime = pm`，`tail` 全 `< dr`。
- `negImpl (term*g)).toList = (⟨dr⟩, -c) :: tail.map neg`，`b := (⟨dr⟩,-c)`，`b.fst.deg = dr`，`b.snd = -c`，
  尾（`negImpl` 不改 fst）全 `< dr`。
- `(a.snd + b.snd).val = (a.snd + (-c)).val = 0`（引理 C，`c.val=a.snd.val, c.prime=pm=a.snd.prime`）。
- 引理 D（`d = dr`）：`mergeAdd (a::rest_r) (b::tail')` 全 `< dr`。即 ALL-DROP。

## 6. 度量下降收尾

`decreasing_by`：目标 `measure r' < measure r`。
`simp only` 展开 measure；`hr : ¬ r.isEmpty` ⇒ RHS `= dr + 1`。
分 `r'.isEmpty`：空 `0 < dr+1`；非空取 `z = r'[0]! ∈ r'.toList`，ALL-DROP 得 `< dr`，`omega`。

## 7. 5 点审核

1. 数学正确性：首项系数 `coeff*lc(g) = a.snd`（乘逆），相减消首项，度数严格降。✓
2. 无跳步：每层模运算展开为 `Nat.mul_mod/mod_mod/mod_eq_of_lt`；分情形 omega。✓
3. Lean 可行：全部 API 已存在（三塔保持引理 + `mergeAdd_lt_all` + `sortedListB_filterMapShift`）。✓
4. 工程：measure 严格降（+1 编码处理 constant÷constant 边界 dr=dg=0：d=0，仍消首项，r' 空或降）；
   `decreasing_by` 所需 bound 全在函数体 `have` 中（不靠 decreasing_by 现推）。✓
5. 边界：r 空/dr<dg 早退不进递归；g 空由 hg_ne 排除；coeff=0 由 NonZeroB 排除。✓

## 8. divmod 调用处

`divmod` 的 decidable 守卫追加 `NonZeroB f ∧ ReducedB f pm ∧ NonZeroB g ∧ ReducedB g pm ∧ g[0]!.fst.deg = dg`，
`pm := g[0]!.snd.prime`，`lc_g_inv := g[0]!.snd.inv`（`h_lc_prime` 需 `inv.prime = pm`——由 inv 定义取同 prime，`rfl` 或小引理）。
其余（`gcdAux/extGcdAux` 是 partial def）不受影响。
