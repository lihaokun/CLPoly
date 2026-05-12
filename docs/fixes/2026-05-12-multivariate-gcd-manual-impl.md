# 修正方案：多变量 GCD 手写实现替代 C++ 翻译路线

> **状态**：草稿，待用户审核。
> 对应 workflow.md §5.1 算法内部错误修正流程（这里是"翻译路线技术受阻"变种）
> 对应分支：feature/formal-proofs
> 前置工作：Phase F-impl-A.1（Math/MvGCD.lean 数学 spec 已就位）

---

## 第一部分：背景与目标

### 上层目标

Phase F-impl 要让 Lean 端 `squarefreefactorize_lex_ir`（已通过 cpp2lean v2 翻译进
Corpus.lean）能 #eval 跑出正确结果，进而支持 B2B 端到端测试。

依赖链：
```
squarefreefactorize_lex_ir (translated)
  ├─ derivative (Model basis)        ← 占位 identity，需真实现
  ├─ polynomial_GCD (Model basis)    ← priority := 0 fallback，需真实现
  ├─ cont (Model basis)              ← 返回常数包装，需真实现多变量
  ├─ leadcoeff (Model basis)         ← priority := 0 fallback，需真实现
  └─ operator/ (Model basis HDiv)    ← identity placeholder，需真实现
```

5 个 basis 原语都需要真实现。

### 两条路线对比

| 路线 | 算法 | Lean 实现 |
|------|------|----------|
| **Y. 翻译**（cpp2lean v2 扩展） | 1:1 C++ Brown's modular GCD | translator emit |
| **X. 手写** | primitive PRS（textbook 算法） | 手写 Lean |

---

## 第二部分：翻译路线（Y）技术评估

### Y 的工作分解

| 阶段 | 内容 | 行数 | 风险 |
|------|------|------|------|
| Y.1 | `survey_ast.py` 支持 multi-instance per function | ~80 Python | 低 |
| Y.2 | `build_pass8_corpus.py` 通用 multi-instance | ~50 Python | 低 |
| Y.3a | **Pass 4 双指针 merge 识别**（数据流分析）| ~150-200 Python | **高** |
| Y.3b | `dense_upoly_zp` 类基础设施（Lean + CLASS_MAP） | ~250 Lean | 中 |
| Y.4 | Pass 5 method 派发 + ctor 处理 | ~80 Python | 中 |
| Y.5 | 各种小补丁 + sorry 清零 | ~100 行 | 中 |
| Y.6 | 67 现有函数回归保护 + 调试 | ~100 行 | **高** |
| **总计** | | **~800-950 行** + **5-8 session** | |

### Y 的关键卡点：双指针并行 merge

C++ `polynomial_GCD` L373-410 段（CRT 增量合并两个有序多项式）：

```cpp
auto Pout_ptr = Pout_.begin();   auto Pout_end = Pout_.end();
auto Pm_ptr   = Pout_mod.begin(); auto Pm_end   = Pout_mod.end();
while (Pout_ptr != Pout_end && Pm_ptr != Pm_end) {
    if (F.comp(Pout_ptr->first, Pm_ptr->first)) {
        tmp_Pout_.push_back({Pout_ptr->first, ...});
        ++Pout_ptr;
    } else {
        if (Pout_ptr->first == Pm_ptr->first) {
            tmp_Pout_.push_back({Pm_ptr->first, ...});
            ++Pout_ptr;
        } else {
            tmp_Pout_.push_back({Pm_ptr->first, ...});
        }
        ++Pm_ptr;
    }
}
// + 两个 tail while 处理剩余
```

**当前 Pass 4 (`iter_recognize`) 处理模式**：
- 单 `for (auto& x : container)` rangefor → `Array.foldl`
- 单 iterator while loop（`while(it != end) {...; ++it;}`）→ rangefor 类似

**双指针 merge 需要**：
- 跟踪两条独立 cursor 状态
- 识别 advance 是条件性的（基于比较）
- emit 双参数递归函数（`mergeAux i j (a b : Array X) : Array Y := ...`）

这是**数据流分析**层级（dataflow analysis），而非 Pass 4 当前的语法模式匹配。需要：
- IR 增加"双 cursor"状态跟踪
- Pass 4 添加 dataflow pattern recognition
- 测试覆盖（merge 在 polynomial add/mul 等多处出现）

虽然有复用价值（识别后所有多项式算术合并都受益），但工程实现的不确定性高。

### Y 的另一个卡点：dense_upoly_zp 类

`__polynomial_GCD` (mod-p 内层) 用 `dense_upoly_zp` 类（dense Zp poly，性能优化）：
- 类型：raw uint64 array + prime
- 方法：constructor、`nmod_mul`、`scalar_mul`、`deg`、`to_upoly`、static `gcd`

Lean 端需新 abbrev + 6+ 方法 + CLASS_MAP entry。约 250 行 Lean basis。

或用 SparsePolyZp（Model 已有）替代 —— 与 C++ 不 1:1，但功能等价。

### Y 成功率估计

| 范围 | 成功率 |
|------|--------|
| Y 完整 1:1 翻译 Brown's | 40-50% |
| Y 部分（cont/leadcoeff 翻译，polynomial_GCD 外层手写） | 70-80% |

---

## 第三部分：手写路线（X）方案

### X 算法选型：primitive PRS

**选 primitive PRS** 而非 subresultant PRS：
- 正确性等价
- 实现更简单（~30% 少代码）
- "中间式膨胀"问题仅对**性能**有影响，对**正确性**没影响
- sqf 调用 GCD 时输入度数不大，膨胀实际影响可控

primitive PRS 算法：
```
gcd(F, G) over R[x] (R = ZZ or ZZ[remaining vars]):
  if F = 0 return G
  if G = 0 return F
  while G ≠ 0:
    R = pseudo-remainder(F, G)
    if R = 0 break
    F = G
    G = primitive_part(R)
  return F
```

**多变量递归**：把多变量视为 univariate over R[remaining vars]，递归至单变量 ZZ。

单变量 ZZ GCD 同样用 primitive PRS。

### X 与 Math/MvGCD.lean spec 的对照

- spec：`gcd_spec` 通过 Mathlib UFD GCDMonoid（noncomputable）
- X 实现：primitive PRS（可执行）
- **数学上等价**：都是 ℤ[x_0, ..., x_n] 的 GCD（associates）
- Refinement proof：在 Phase F-proof 阶段证 `gcdImpl = gcd_spec`（模 associates）

注：与 C++ 的 Brown's modular GCD 不是 1:1 算法，但与 Mathlib spec 数学等价。

### X 工作分解

| 阶段 | 内容 | 行数 |
|------|------|------|
| X.1 | `MvPolyZZ.derivativeMv` 真实现（C++ 1:1） | ~30 |
| X.2 | `MvPolyZZ.leadcoeffMv` 真实现（C++ 1:1） | ~30 |
| X.3 | 辅助：multivariate 主变量 deg / coef-by-deg 分组 | ~50 |
| X.4 | `MvPolyZZ.contMv` 真实现（用 polynomial_GCD 递归）| ~80 |
| X.5 | `MvPolyZZ.polynomialGCD`（primitive PRS） | ~180 |
| X.6 | `SparsePolyZZ.gcd`（单变量 ZZ 的 PRS GCD） | ~80 |
| X.7 | `HDiv MvPolyZZ MvPolyZZ MvPolyZZ`（exact division，leading-mono reduction） | ~50 |
| X.8 | 测试：sqf 端到端 + 各原语单元测试 | ~80 |
| **总计** | | **~580 行** |

### X 成功率

**95%+** —— 经典算法，无翻译器风险，无 dataflow 分析。

时间：**2-3 session**（含测试和调试）。

---

## 第四部分：路线决策与权衡

| 维度 | Y（翻译）| X（手写）|
|------|---------|---------|
| 工作量 | 800-950 行 + 5-8 session | 580 行 + 2-3 session |
| 成功率 | 40-50%（完整）/ 70-80%（部分）| 95%+ |
| 1:1 C++ Brown's | ✓（如成功）| ✗（不同算法）|
| Lean 真实现 + B2B 真 cross-check | ✓ | ✓ |
| 复用价值 | 高（双指针 merge / dense_zp 类，未来其他翻译受益）| 低 |
| Refinement 证明对照 | ✓（Mathlib spec）| ✓（同 spec，但算法不同）|
| 违反项目"1:1 C++"原则 | ✗ | ✓（需此文档说明）|

### 关于 "1:1 C++" 原则的违反

`proof/CLAUDE.md §L2 算法模型原则`：
> L2 模型必须 1:1 对应 C++ 算法逻辑

X 路线选 primitive PRS（不是 C++ Brown's modular GCD），**严格说违反**。

理由（本文档存在的意义）：
1. Y 路线技术受阻：双指针 merge + dense_upoly_zp 类是高风险扩展
2. X 路线快速实现：B2B 端到端测试需求紧迫
3. **数学等价**：都正确求 GCD，与 Mathlib spec 一致
4. **未来可扩**：Y 路线可作为后续 Phase F-impl-v2，把 X 替换为翻译产物

### Refinement proof 影响

- Math/MvGCD.lean 的 spec（A.1 已就位）保持不动
- L1 实现（X 路线手写）需独立 refinement proof：`MvPolyZZ.gcdImpl ≈ gcd_spec`（associates）
- 若未来切到 Y 翻译，原 refinement proof 需重做（翻译产物算法不同）

---

## 第五部分：测试计划

### 单元测试

| 测试 | 内容 |
|------|------|
| `test_derivative_mv` | x², x²y, -2x³+xy 等求导 |
| `test_leadcoeff_mv` | 主变量首项系数（剩余变量中的多项式） |
| `test_cont_mv` | 整数 cont（单变量）+ 真多变量 cont |
| `test_gcd_mv_uni` | univariate-as-mvpoly GCD：gcd(x²-1, x-1) = x-1 |
| `test_gcd_mv_truly` | 真多变量：gcd((x+y)², x+y) = x+y |
| `test_hdiv_mv` | exact division: (x²-1)/(x-1) = x+1 |

### 端到端：sqf 测试

| 输入 | 期望 |
|------|------|
| `x` (univariate) | `[(1,1), (x,1)]` |
| `x²` | `[(1,1), (x,2)]`（重数 2）|
| `(x-1)²·x` | `[(1,1), (x,1), (x-1,2)]` |
| `-6(x-1)(x-2)` | `[(-6,1), (x²-3x+2,1)]`（issue #14 兼容）|

### 全量回归

- `lake build` 3071+ 全绿
- B2B 35/36 PASS 不退（加 sqf 向量）

---

## 第六部分：执行步骤（待用户批准）

1. **[等批准]** 用户确认 X 路线
2. X.1 (derivative) + 测试 → [报告]
3. X.2 (leadcoeff) + 测试 → [报告]
4. X.3-X.4 (主变量辅助 + contMv) + 测试 → [报告]
5. X.6 (SparsePolyZZ.gcd) + 测试 → [报告]
6. X.5 (polynomialGCD multivar) + 测试 → [报告]
7. X.7 (HDiv) + 测试 → [报告]
8. X.8 端到端 sqf 测试 + B2B 向量 + 全量回归
9. devlog + commit + push

---

## 第七部分：风险与未决问题

1. **Q**: primitive PRS 对深嵌套多变量（如 4+ 变量）会有中间式膨胀，性能差。会不会让 B2B 跑超时？
   **A**: sqf 的实际使用主要是单变量 + 少量浅多变量，影响有限。若实际 B2B 超时，再优化。

2. **Q**: 多变量主变量选哪个？C++ 用 `get_first_var(F)`。Lean Model 的 lex order 是 var ID 升序还是降序？
   **A**: 需在实施时核对（看 Monomial.lt 定义）。

3. **Q**: 单变量 ZZ GCD（SparsePolyZZ.gcd）是否要也用 modular？还是直接 PRS？
   **A**: PRS。Model 已有 SparsePolyZp.gcd（mod-p），但单变量 ZZ 直接 PRS 简单且对正确性足够。

4. **Q**: 与 Phase F-impl-A.1（Math/MvGCD.lean spec）的关系如何描述？
   **A**: spec 是数学等价类，impl 是具体算法。Refinement proof 在 Phase F-proof 中证它们相等（模 associates）。
