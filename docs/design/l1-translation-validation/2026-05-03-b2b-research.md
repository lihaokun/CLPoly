# B2B 语义测试 — 调研报告

**日期**：2026-05-03
**前置**：cpp2lean v2 翻译器 lake build 100% 通过（commit `cbfd7aa`），但 Lean Model 中大量 stub 是 `default` 占位
**目标**：评估 B2B（back-to-back）测试可行性，给出实施起点

## 1. 现状摸底

### 1.1 翻译器状态

- 346 函数全部 lake build 通过（类型检查✓，语义✗）
- Lean Model（`proof/lean/CLPoly/Model.lean`）共 134 个 def，~80 个是 stub
- stub 类型：
  - `:= default` 直接占位（10+）
  - `(_x : α) : β := default`（多态 stub）
  - `f a b := default` / `_a := target`（identity 占位）

### 1.2 现有 B2B 设计

`docs/design/l1-translation-validation/back2back-design.md` 已存在（164 行），但：
- 只覆盖 `polynomial_factorize_zp.hh` 13 个 Zp 函数
- 当前 corpus 是 346 函数（含 univar/multivar/Wang/MTSHL）
- 设计中 "翻译后 IR 文件" 词汇与当前 v2 实际不一致（v2 emit 单个 `Generated/Corpus.lean`）

### 1.3 #eval 可行性（实测）

| 测试 | 结果 |
|------|------|
| `__make_zp_ir 7 13` | ✅ `{val := 7, prime := 13}` |
| `__upoly_make_monic_ir poly` | ✅ 正常返回 |
| `__squarefree_Zp_ir poly`（依赖 stub `polynomial_GCD := default`） | ⚠️ 返回 `#[]`（语义错，但 #eval 不挂） |

结论：
- **#eval 可用**（partial def 不挂）
- **stub 让结果错** —— 必须先有真实 Lean 实现的原语

## 2. B2B 测试关键问题

### 2.1 stub 替换 vs. opaque + spec 等价

两种思路：

**A. stub 替换为可计算 Lean 实现**（设计文档当前路径）：
- `polynomial_GCD` → 用 Mathlib `EuclideanDomain.gcd` 或手写欧几里得
- `derivative` → 手写求导
- 优点：能 `#eval`，直接对比 C++ 输出
- 缺点：~30+ 原语需手写实现；某些（`squarefreefactorize` 多变量）实现复杂

**B. opaque + spec 等价**（理论路线）：
- 把原语标 `opaque`，证明翻译保留语义而非数值
- 优点：不写实现；和 L2 证明思路一致
- 缺点：不能 `#eval`，需要纯证明；本阶段（验翻译忠实性）不适用

**结论**：走 A 路线，但分批——先实现核心原语（够测试 Zp/单变量），多变量原语（squarefreefactorize 等）后续。

### 2.2 测试向量来源

| 来源 | 优劣 |
|------|------|
| 手写 JSON | 可控但覆盖窄 |
| 随机生成 | 覆盖广但可能踩 require fail（边界） |
| C++ test suite 录制 | 真实场景，但需 C++ 端改造 |

**建议**：起步阶段手写 + 现有 C++ test cases 提取（小 deg + 小 prime 用例）。

### 2.3 函数选择（覆盖优先级）

不是 346 函数都同等重要。按依赖深度分层：

| 层 | 函数（示例） | stub 依赖 |
|----|------------|----------|
| **L0 原子** | `__make_zp`, `Zp.ofInt`, `get_deg`, `Array.head!` | 无 stub（直接计算）|
| **L1 单变量基础** | `__upoly_make_monic`, `__upoly_subtract_x`, `__upoly_mod` | 依赖 L0 |
| **L2 单变量算法** | `__squarefree_Zp`, `__ddf_Zp`, `__edf_Zp` | 依赖 `polynomial_GCD`, `derivative`（需实现）|
| **L3 单变量顶层** | `__factor_Zp`, `factorize_upoly_ir` | 依赖 L2 |
| **L4 多变量** | `__factor_multivar`, `__wang_core`, MTSHL 系列 | 依赖 `squarefreefactorize`, `assign`, `poly_convert` 等大量 stub |

**起点**：从 L0/L1 开始（无 stub 依赖，直接可测），逐层上爬。

### 2.4 比较粒度

| 粒度 | 说明 |
|------|------|
| 数值精确 | `(0 : Int)` vs `(0 : UInt64)` 算 fail？ |
| 语义等价 | 多项式 `2x + 1` 和 `1 + 2x` 算同（对应顺序无关）？ |
| 容器忽略 | Factorization 的 factors 顺序？ |

**建议**：先精确比较，遇到顺序无关情况单独 normalize。

## 3. 实施起点（建议）

按工作流，分 3 步：

### 步骤 1：实现核心原语（~1 天）

实现 8-10 个 L0/L1 原语的真实 Lean 算法：

| 原语 | 实现 |
|------|------|
| `polynomial_GCD` | 单变量欧几里得 GCD（10-20 行）|
| `derivative` | 多项式求导（~10 行）|
| `__upoly_mod` | 多项式取模 |
| `Zp.invert` | 扩展欧几里得求模逆 |
| `pair_vec_div` | 多项式长除法 |
| `polynomial_mod` | 系数 mod p |
| `cont` / `pp` | content / primitive part |
| `get_first_deg` | 首项度数 |

完成后能 `#eval` L1+L2 函数返回有意义的结果。

### 步骤 2：B2B 驱动框架（~1 天）

- C++ 端：`test/test_b2b.cc` 调用 13 个 Zp 函数（先按现有 design），输出 JSON
- Lean 端：生成 `eval_b2b.lean` `#eval` 同一组测试向量
- 比较脚本：`tests/b2b_compare.py`
- 测试向量：手写 50-100 个（覆盖 L0-L3）

### 步骤 3：跑通 + 修翻译错误（不定）

预计前几次跑会发现：
- 翻译器 bug（参数顺序错、cast 漏等）→ 修 Pass 1-8
- Lean Model stub 实现 bug → 修
- 测试向量 require 不满足 → 改

## 4. 不确定项

1. **`partial def` 死循环**：~Wang/MTSHL 函数实际可能不终止（C++ 端 fuel-bound，Lean 无 fuel）。需要 timeout + 跳过。
2. **多变量原语**（squarefreefactorize 等）：实现复杂，可能阶段 2 才做。
3. **数值表示差异**：C++ ZZ 和 Lean Int / GMP；C++ Zp 用 UInt64，Lean Zp.val 用 UInt64。理论一致但要验。
4. **结果序列化**：JSON 格式 + Lean 端打印格式对齐（如 `Zp.mk` vs `{val := X, prime := Y}`）。

## 5. 推荐起点

**今天可以启动的工作**：

1. 把 `back2back-design.md` 更新到 v2 翻译器现状（13 函数 → 346 函数分层）
2. 实现核心 8-10 个原语（步骤 1）
3. 用 4-5 个 L0/L1 函数搭起 B2B 最小闭环（端到端跑通一个 PASS）

**不推荐**：
- 一开始就攻 multivar（依赖太多 stub，调试地狱）
- 一次性写所有测试向量（先跑通最小闭环再扩）

## 6. 决策点

请确认以下方向：

1. **stub 策略 A（实现可计算 Lean 原语）vs B（opaque + spec）**：建议 A
2. **起步范围**：13 个 Zp 函数 vs 346 全部 vs 手选 L0-L1（~30 个）：建议 L0-L1 起步
3. **测试向量**：手写起步 vs 直接随机：建议手写 50-100 个
4. **C++ 端修改**：写新 `test_b2b.cc` vs 复用 `test_crosscheck_*.cc`：建议新写（专注翻译验证）
5. **是否更新 back2back-design.md**：建议先更新（设计三阶段第 1 步）

确认后进入"架构"阶段，写详细方案。
