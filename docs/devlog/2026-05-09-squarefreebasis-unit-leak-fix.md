# 2026-05-09: squarefreebasis unit leak 修复（issue #14）

## 做了什么

修复 GitHub issue #14：`squarefreebasis` 在输出 basis 中混入常数 `-1`，且数量随输入规模增长。

**根因**：`squarefreefactorize` 没把 factor 首项 sign 提到 cont，导致同一个数学多项式（associates）在 `squarefreefactorize` 输出和 `polynomial_GCD` 输出之间字面表示不齐 → `squarefreebasis` 主循环用 `F1 != lst[j]` 字面比较时误判为"应拆分" → 拆分时 `lst[j]/F1` 是 unit 残留进 basis。

**修复**：在 `clpoly/polynomial_gcd.hh:109-111` 的 `squarefreefactorize` 中加 4 行 sign 归一化：当 `F_ = F/cont(F)` 的首项系数为负时，把 sign 从 F_ 转移到 f_cont（保持 cont * F_ = F 不变）。这样 factor 永远首项正，与 `polynomial_GCD` 的内部归一化对齐。

## 为什么做

- 用户报告的 user-visible bug，影响 charset 的 squarefree 基底结果
- 根因不只是 squarefreebasis 局部问题，而是 sqf 的 sign convention 与业界（SymPy/Maple/FLINT）不一致
- 在 sqf 上修是治根方案，sqfb 自然就对（关联代码完全不用动）

## 关键决策及其理由

### 决策 1：选择 Plan B（上游 sqf 修复）而非 Plan A（局部 sqfb 修复）

- **A** 是 5-10 行局部判定补丁
- **B** 是 4-5 行上游 sign-normalize
- B 完成后 A 是 dead code（不会触发的判定）
- B 一次到位，避免在 codebase 留"已知冗余的局部补丁"
- B 让 CLPoly 的 sqf 输出与 SymPy/Maple/FLINT 完全一致

### 决策 2：cont sign convention 选"跟随 leading"（不动 cont 函数）

- 不改 `cont` 函数本身（保持现有 GMP gcd 行为：多项情形返回正 gcd）
- 仅在 sqf 内部处理 sign：`F_ = F/cont` 后若首项负则双方取负
- → 不影响其他调用 `cont` 的代码

### 决策 3：调用点扫描（避免破坏链路）

实测 7 处 `squarefreefactorize` 调用点 + 2 处测试断言，全部对 factor sign 不敏感：
- `charset.hh`: `G * factor` 乘积语义
- `polynomial_factorize*.hh`: `is_number` 检查或递归 factor
- `test_polynomial_gcd.cc`: `divides` + `is_squarefree` 断言（不依赖 sign）

→ 修复无需改任何调用方。

## 遇到的问题与解决方式

### 问题 1：第一次直接动代码违反 §5.1 流程

第一次尝试时跳过了"先写修正方案 → 用户确认 → 再实施"流程，直接 patch + 测试，被用户纠正。

**解决**：
1. revert 全部代码改动
2. 严格按 workflow.md §5.1 流程：先写 `docs/fixes/` 修正方案文档（含 SymPy/Maple/FLINT 对照 + 调用点扫描）
3. 等用户审核 + 确认方案后再实施

### 问题 2：第一版修正方案文档过于复杂

第一版列了 Plan A/B/C 三方案对比 + 模糊推荐 C（先 A 后 B），用户指出："先 B 后 A 不更合适？"

**解决**：重写文档为纯 Plan B 路线，把 A/C 全部删除，理由也改写为"A 在 B 完成后是 dead code"。简化决策路径。

### 问题 3：错误评估 sympy 一致性

第一次声称"CLPoly 与 SymPy 不一致"，后实测 SymPy 输出后发现仅在「多项 + 首项负」一个具体场景下不一致，其他完全一致。这让 Plan B 的工作量从估计的 30-50 行降到实际的 4 行。

**教训**：参考实现对照必须**实测 ground truth**，不能靠记忆/文档猜测。

## 涉及的文件

### 主修改
- `clpoly/polynomial_gcd.hh`：+8 行（4 行修复 + 4 行注释），位于 `squarefreefactorize` L109-111 周围

### 新增文件
- `docs/fixes/2026-05-09-squarefreebasis-unit-leak.md`：修正方案文档（§5.1 必需）
- `test/test_squarefreebasis_bugfix.cc`：回归测试（10 用例覆盖 issue 复现 + canonical 断言 + sanity）
- `docs/devlog/2026-05-09-squarefreebasis-unit-leak-fix.md`：本日志

### 流程改进
- `test/run_all_tests.sh`：把 `test_factorize_bugfix` 和 `test_squarefreebasis_bugfix` 加入 TESTS 数组（之前 bugfix 测试都没进 runner，是历史疏漏）

## 度量

- 耗时：~3 小时
  - 第一次直接动代码：~30 分钟（被 revert）
  - 修正方案 v1（Plan A/B/C 对比）：~30 分钟
  - 修正方案 v2（Plan B 路线 + sympy 实测 + 调用点扫描）：~45 分钟
  - 实施 + 测试：~1 小时
  - devlog + commit：~15 分钟
- 迭代：~5 轮（含 v1 文档被指正、SymPy 实测后认知修正等）
- C++ 修改行数：12 行（含注释）
- 新增测试代码：~140 行（10 用例）
- 修正方案文档：~200 行
- 放弃的方案：
  - **Plan A（局部 sqfb 内 quotient is_number 检查）**：被 B 吞掉，不做
  - **Plan C（先 A 后 B 双管）**：A 部分被 B 吞掉，等于浪费一轮 commit
  - **Plan B 30-50 行重写 cont sign convention**：实测后发现只需 4 行 sqf 局部 sign 归一化即可达到目标，不动 cont
