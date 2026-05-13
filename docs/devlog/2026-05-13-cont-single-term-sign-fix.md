# 2026-05-13: cont 单项情形 sign 不一致修复（issue #17）

## 做了什么

修复 GitHub issue #17：CLPoly 的 `cont()` 函数对单项与多项多项式返回不一致的 sign，违反 Maple/SymPy/FLINT 通用"content always non-negative"约定。

3 处改动：
- `clpoly/polynomial_gcd.hh:1067`（univariate ZZ cont 特化版）
- `clpoly/upolynomial.hh:220`（template Tc cont 版本）
- `clpoly/polynomial_gcd.hh:581+`（多变量 cont 函数）

## 为什么做

发现于讨论 Phase F-impl 时审视 cont sign 行为：用户问 "cont(-x) = -1 是 CLPoly 还是 master 的行为？" 实测确认 CLPoly C++ 单项 cont 直接返回 leading coef（带符号），多项 cont 通过 GMP `mpz_gcd` 自动折正。不一致。

进一步对照 Maple/SymPy/FLINT 确认这是 CLPoly 独有 quirk，违反工业约定。

issue #14 (squarefreebasis unit leak) 的根因之一也是这个：PR #16 在 squarefreefactorize 加 sign-normalize 兜底修了表象，但 cont 不一致根因未修。

## 关键决策及其理由

### 决策 1: scope 包含 univariate + multivariate cont（不分两 PR）

初版 doc 写"多变量留 follow-up issue #17b"，被用户指出是 artificial split。修订为同 PR 修两处，理由：
- 是同一个 bug（cont sign 约定）
- 修复方式独立（univariate: 改 init c=0；multivariate: 末尾加 normalize 3 行）
- 风险隔离意义不大（pp 等下游影响同时存在）

### 决策 2: 半形式化证明

用户 review 时指出文档缺正确性论证。补：
- Spec 数学定义（含 CLPoly "cont(0)=1" 历史约定与 Maple "cont(0)=0" 差异说明）
- Lemma 1（loop invariant，归纳证明）
- Main Theorem（实现 = Spec，分 case）
- 等价性 Claim（修后多项情形与原实现等价，证明）
- 边缘情形枚举表（5 行覆盖空/单项正/单项负/单项零/多项）

### 决策 3: 修复模式选择

univariate cont 候选：
- A. 初始 `c = abs(first coef)`，然后正常 gcd 折叠（保留原循环结构）
- B. 初始 `c = 0`，range-for 折叠所有项（更简洁，迭代统一）

选 B。两种等价（首步 `gcd(0, c_0) = |c_0|`），B 代码更短更清晰。

multivariate cont：
- 在 return 前加 sign normalize 3 行（cover single-group 和 multi-group 两条路径，对已归正情形是 no-op）
- 比逐处修补"cont = std::move(tmp)"分支更鲁棒

## 遇到的问题与解决方式

### 问题 1: 第一版 doc scope 不够

初版只修 univariate，留多变量为"未来 issue #17b"。被用户质疑是 artificial scope split。修订为同 PR 修两处。

### 问题 2: 文档缺正确性证明

用户问"你进行半形式化证明了吗？"，doc 之前只有非正式"`gcd(0,x)=|x|` 所以对"。补正式 Lemma + Theorem。

### 问题 3: 测试 pp 调用时模板序不匹配

`polynomial_<ZZ>` 默认 grlex 序，pp 定义在 `polynomial_<ZZ, lex>` 上。测试 pp invariant 编译失败。

解决：去掉 pp 调用，改测 cont 非负 + 用 lex 序明确做多变量 cont 测试。

### 问题 4: 现有测试 expected 值改变

预期是有的，但实际所有现有测试（含 test_polynomial_gcd, test_upolynomial, test_multivar_helpers, test_factorize_trace）全过。因为它们的断言都是"乘积等于"或"非空"等 sign-不敏感谓词。

## 涉及的文件

### C++ 改动
- `clpoly/polynomial_gcd.hh:1067` ZZ 特化版 cont（重写 - init c=0）
- `clpoly/polynomial_gcd.hh:631+` 多变量 cont（末尾加 sign normalize 4 行）
- `clpoly/upolynomial.hh:220` template Tc cont（重写）

### 新增文件
- `docs/fixes/2026-05-12-cont-single-term-sign.md` 修正方案文档（416 行）
- `test/test_cont_sign_consistency.cc` 回归测试（16 用例）
- `docs/devlog/2026-05-13-cont-single-term-sign-fix.md` 本日志

### 流程
- `test/run_all_tests.sh` 加 `test_cont_sign_consistency` 到 runner

## 度量

- 耗时：~4 小时
  - 实测 + 调研 C++/Maple/SymPy/FLINT 对照：~30 min
  - 修正方案 doc 第一版（含 multivariate 误判为 follow-up）：~30 min
  - 用户 review 反馈 → doc 修订（含 半形式化证明）：~1 hour
  - 二轮 self-audit：~30 min
  - 实施 + 测试 + 调试：~1 hour
  - devlog + commit：~15 min
- 迭代：~3 轮（含 doc scope 修订、证明补充、测试调试）
- C++ 改动行数：~25 行（含注释）
- 新增测试代码：~140 行
- 修正方案文档：416 行（含半形式化证明）
- 放弃的方案：
  - 多变量 cont 留 follow-up（被用户指出 artificial split）
  - 替代修复模式"abs(first) + 保留循环"（选了更简洁的"c=0 + range-for"）
  - pp 联动测试（模板序不匹配，简化掉）
