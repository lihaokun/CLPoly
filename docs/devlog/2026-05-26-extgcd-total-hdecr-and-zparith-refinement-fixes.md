# extGcdAux totality & ZpArith refinement fix

Date: 2026-05-26

## 做了什么

**Model.lean:**
- `extGcdAux` 改为 total `def`，补全 `h_decr` 终止性证明
  - `hq` calc 块简化（去掉嵌套 calc）
  - `exact_mod_cast` → `Int.ofNat_lt.mp`（纯 Lean 4 无 `mod_cast`）
  - 补回遗漏的 `exact h_natAbs_lt`（之前 edit 时丢掉了这行）

**ZpArith.lean (`__upoly_make_monic_ir_refines`):**
1. 空数组分支：`simp` 无法处理 `SparsePolyZp.front! #[]`（= `default`），改用 `Array.isEmpty_iff` → `subst` + `by_cases` 分支
2. `h_one` 分支：`h_modInv_eq_one` 直接用 `modInv_one` 构造（而不是先 `rw [h_val]` 再调用），避免 `simp` 重写顺序问题
3. `¬hone_val` 分支：`hne_beq` 是 `==`（Bool），目标条件是 `=`（Prop），需要 `h_ne_prop` 用 `hne_val` + `native_decide` 转型

## 涉及文件

- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Refinement/ZpArith.lean`
