# nl-proof 草稿：`__squarefree_Zp_ir_refines` Branch B（导数≠0）

目标（`Refinement/SquarefreeZp.lean:2371` 的 admit）：在强归纳的 `ih` 下，导数≠0 分支证
```
toPolyList (Generated.__squarefree_Zp_ir g) p = sqfZp (SparsePolyZp.toPoly p g)
```
其中 `F := SparsePolyZp.toPoly p g`，`h_deriv0 : ¬ (derivative g = empty)`，`h_deg0_g : ¬ F.natDegree = 0`。

## 1. 两侧结构对照

**L2 `sqfZp F`（else-else 分支）**
```
c  := normalize (gcd F F')            -- F' = derivative F
w  := normalize (F /ₘ c)
(yun_result, c_rem) := yunLoop w c 1 [] hc_ne
if 0 < c_rem.natDegree then
  yun_result ++ (sqfZp (contract p c_rem)).map (fun (h,e) => (h, e*p))
else
  yun_result
```

**C++ `__squarefree_Zp_ir_def g`（else 分支，Corpus.lean:224-247）**
```
c_1 := polynomial_GCD g (derivative g)                 -- 未 normalize 的 gcd
w_3 := normalization (pair_vec_div _ g c_1 (comp g))   -- normalize (g /ₘ c_1)
(_, c_2, result_4) := _loop___squarefree_Zp_1_ir_def 1 w_3 c_1 #[]   -- C++ yun 循环
if (¬isEmpty c_2 ∧ get_deg c_2 > 0) then               -- 0 < c_rem.natDegree
  g_3 := __extract_pth_root_ir_def c_2                 -- contract p c_rem
  g_4 := (__upoly_make_monic_ir_def g_3).snd           -- normalize (contract p c_rem)
  sub_3 := __squarefree_Zp_ir_def g_4                  -- 递归
  result_7 := (_loop___squarefree_Zp_2_ir_def 0 sub_3 result_4 p_1).2.2  -- result_4 ++ sub_3.map(×p)
  result_7
else
  result_4
```

**C++ yun 循环 `_loop___squarefree_Zp_1_ir_def i w c result`（Corpus.lean:182-195）**
```
if (¬isEmpty w ∧ get_deg w > 0) then          -- ¬ w.natDegree = 0
  y := polynomial_GCD w c                       -- 未 normalize
  z := normalization (pair_vec_div _ w y (comp w))   -- normalize (w /ₘ y)
  if (¬isEmpty z ∧ get_deg z > 0) then          -- 0 < z.natDegree
    result' := push (z, i) 到 result             -- acc ++ [(z,i)]
    (走 bb_14)
  else result' := result
  bb_14: c' := normalization (pair_vec_div _ c y (comp c))  -- normalize (c /ₘ y)
         recurse i+1, w:=y, c:=c', result'
else (0, c, result)                             -- (acc, c)
```

**关键差异（根因级，必须在证明中弥合）**：
- C++ 的 `y = polynomial_GCD w c` **未 normalize**；L2 `y = normalize(gcd w c)`。
- 故 C++ 每层递归携带的 `w := y`（未 normalize）与 L2 携带的 `w := normalize(gcd)` **只相伴（associated），不相等**。
- 但两侧输出的 `z`、`c'` 都外套 `normalize`，在 w,c 相伴时相等。
- ⇒ 需要一条**相伴不变量**贯穿循环归纳。

## 2. 前置桥接引理（按依赖顺序，含现状）

| 编号 | 引理 | 陈述（toPoly 语义） | 现状 | 策略 |
|------|------|--------------------|------|------|
| P1 | `divmod_deg_decrease`(445) | `r' = r - term*g ⇒ deg r' < dr` | admit | 首项抵消：`term` 使 `r'` 的 dr 次系数为 0；用 `SparsePolyZp` 加/乘的 toPoly 语义 + `natDegree` 单调 |
| P2 | `divmod'_wf_deg_lt`(1081) | `deg (toPoly (divmod'_wf f g).snd) < deg (toPoly g)`（g 非空） | admit | 由 P1 + `divmodAux'_wf` 展开；余项每步 dr 递减、终止时 dr<dg |
| P3 | `divmodAux'_wf_eq`(692) | WF 版 = partial def 版 | admit | **难点**：partial def 不可归纳。方案：改为让 `pair_vec_div`/`polynomial_GCD` 的 toPoly 语义**直接**用 WF 模型证（P4/P5），绕开 partial-def 等价——即证 `toPoly (C++ divmod) = toPoly f - toPoly g * q` 形式的 divmod 恒等式（已有 `divmod'_wf_identity`@906），不依赖 P3 |
| P4 | `gcdAux'_wf` 终止(673) | measure `deg(toPoly g)` 递减 | admit | 由 P2 直接得 |
| P5 | `polynomial_GCD_toPoly` | `Associated (toPoly (polynomial_GCD f g)) (gcd (toPoly f) (toPoly g))` | 新建 | `polynomial_GCD = gcdAux'_wf`（HasPolyGCD 实例）→ 对 `deg(toPoly g)` 归纳，用 `divmod'_wf_identity`（余项 = f - g*q）得 `gcd f g = gcd g (f mod g)`，对应 Euclid |
| P6 | `pair_vec_div_toPoly` | `toPoly (pair_vec_div _ f g _) = toPoly f /ₘ toPoly g`（g monic 时精确） | 新建 | 用 `divmod'_wf_identity` + `div_modByMonic_unique` |
| P7 | `normalization_toPoly` | `toPoly (SparsePolyZp.normalization f) = normalize (toPoly f)` | 新建 | normalization = 除以首项系数 ⇒ `C (lc⁻¹) * f`，`normalize` 亦然（域上 `normalize f = C (leadingCoeff f)⁻¹ * f`） |

## 3. 主循环引理 B1（yun 循环精化）

**陈述**（相伴不变量归纳）：设 C++ 状态 `(i, W, C, R)` 与 L2 状态 `(w, c, i, acc)` 满足
- `Associated (toPoly p W) w`，`Associated (toPoly p C) c`（且都 WellFormed/AllReduced/val≠0）
- `toPolyList R p = acc`（列表逐项相等）
- 相同 `i`
则
```
(令 (_, C_out, R_out) := _loop___squarefree_Zp_1_ir_def i W C R
    (acc_out, c_out) := yunLoop w c i acc hc)
Associated (toPoly p C_out) c_out ∧ toPolyList R_out p = acc_out
```

**证明**：对 `w.natDegree`（或 `W.natDegree`）强归纳。
- **base**：`w.natDegree = 0`。C++ 侧 `¬isEmpty w ∧ deg>0` 为假（因 `deg W = deg w = 0`，需 `Associated ⇒ natDegree 相等` + `get_deg` 对应 natDegree），返回 `(C, R)`；L2 返回 `(acc, c)`。由不变量成立。
- **step**：两侧都进入递归。
  - `y_C := polynomial_GCD W C`，`y_L := normalize(gcd w c)`。由 P5 + 不变量相伴 + `sqfZp_smul` 用过的 gcd 相伴技巧 ⇒ `Associated (toPoly y_C) y_L`。
  - `z_C := normalization(W /ₘ y_C)`，`z_L := normalize(w /ₘ y_L)`。由 P6+P7 + 相伴 ⇒ `toPoly z_C = z_L`（外套 normalize 消去相伴单位，**精确相等**）。故 push 条件 `0<deg z` 两侧一致，push 的 `(z,i)` 相等。
  - `c'_C := normalization(C /ₘ y_C)`，`c'_L := normalize(c /ₘ y_L)`：同理 `toPoly c'_C = c'_L`（精确相等 ⇒ 也相伴）。
  - 新 `W := y_C`，`w := y_L`：`Associated (toPoly y_C) y_L` 成立（上一步）。
  - 度数下降：`deg y < deg w`（gcd 整除 w，且 w deg>0）——复用 L2 `yunLoop` 的 `decreasing_by` 论证（已在 Algorithm 里证过，可搬）。
  - 应用 ih。

**注**：C++ 的 `w := y_C` 与 L2 `w := y_L` 只相伴不相等，是本引理必须用相伴不变量（而非等式）的原因。`c'`、`z` 因外套 normalize 而精确相等。

## 4. 组装 B2

1. `rw [if_neg h_deg0_g]` 把 `__squarefree_Zp_ir_safe`（或 `__squarefree_Zp_ir`）展开到 C++ else 分支。
2. C++ 初始：`c_1 = polynomial_GCD g g'`，`w_3 = normalize(g /ₘ c_1)`。
   L2 初始：`c = normalize(gcd F F')`，`w = normalize(F /ₘ c)`。
   - `Associated (toPoly c_1) c`：P5（`gcd g g'` 相伴 `normalize(gcd F F')`，注意 `toPoly g = F`、`derivative` 桥接用已有 `derivative_toPoly_eq`）。
   - `toPoly w_3 = w`：`normalize(g /ₘ c_1)` vs `normalize(F /ₘ normalize(gcd F F'))`——除数 `c_1` 与 `normalize(gcd)` 相伴，外套 normalize 消单位 ⇒ 相等（同 §3 技巧）。
3. 应用 B1（初始不变量：`R=#[]↔acc=[]`，i=1）⇒ `Associated (toPoly c_2) c_rem ∧ toPolyList result_4 p = yun_result`。
4. `if 0 < get_deg c_2`：由 `Associated ⇒ natDegree 相等` 对齐 L2 的 `if 0<c_rem.natDegree`。
   - **then**：
     - `g_3 = extract_pth_root c_2` ↔ `contract p c_rem`：用 `extract_pth_root_toPoly_eq`（已证），但需 `toPoly c_2` 与 `c_rem` 关系——只相伴！`contract` 对相伴的输入…需 `contract p` 尊重相伴？`contract` 非乘法同态，`Associated ⇒ contract 相伴`需单独证（或用 `c_2` 与 `c_rem` 都 normalize 后相等：实际 C++ `c_2` 是循环输出的 normalize 结果，L2 `c_rem` 亦然 ⇒ 可能精确相等，需核查 B1 输出是否精确相等而非仅相伴）。**待定：核查 c_rem 是否精确相等**。
     - `g_4 = make_monic g_3` ↔ `normalize`，递归 `sub_3 = __squarefree_Zp_ir g_4`：用 ih（度数下降 `deg(contract c_rem) < deg F`，同 Branch A/L2 decreasing_by）。
     - `result_7 = result_4 ++ sub_3.map(×p)`：用 `h_loop_lemma` 类比（`_loop___squarefree_Zp_2` 与 Branch A 的 `_loop___squarefree_Zp_0` 同构）+ **`h_toPolyList_specific` 溢出界（2354，乙类）**——此处再次出现，需乙类不变量。
     - `sqfZp_smul`（已证）处理 make_monic 的单位。
   - **else**：直接 `result_4 ↔ yun_result`。

## 5. 未决 / 与其它 admit 的耦合

- **乙类 2354 溢出**：Branch B 的 then 分支同样需要"输出乘数 × p < 2⁶⁴"。⇒ 建议在 B1/B2 里一并维护"输出乘数 ≤ natDegree"不变量（见 `docs/fixes/…overflow-precondition.md`）。
- **P3 (`divmodAux'_wf_eq`)**：若能用 `divmod'_wf_identity` 直接给出 P5/P6 的 toPoly 语义，则可**绕开** partial-def 等价这个最难 admit。需核实 `pair_vec_div`（HasPolyDivmod 实例）的 toPoly 是否已有直接引理，或必须过 `divmodAux'_wf_eq`。
- **c_rem 精确相等 vs 相伴**：§4 then 分支依赖此，需在 B1 中确认 `c_out` 是否精确相等（C++ 与 L2 都对 c' 外套 normalize，倾向精确相等；若是则 §4 简化）。

## 6. 建议实现顺序

P2←P1 → P4 → P5,P6,P7（其中设法绕开 P3）→ B1（yun 循环，最大）→ B2（组装，含调用 ih + 乙类不变量）。
估计：~500–800 行 Lean。P5/B1 为主要瓶颈。
