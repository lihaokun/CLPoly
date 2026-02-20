# Gröbner 基调研报告

> 目标：为 CLPoly 添加 Gröbner 基计算支持。
> 本报告梳理 CLPoly 已有基础设施，调研外部主流实现，评估算法方案。

---

## 1. CLPoly 已有基础设施

### 1.1 可直接复用的组件

| 组件 | 位置 | 说明 |
|------|------|------|
| 单项式序 `lex`、`grlex` | `monomial_order.hh:70,95` | Gröbner 基依赖单项式序 |
| 模板化序支持 | `basic_monomial<compare>` | 同一实现可支持任意序 |
| 多项式算术 `+`、`-`、`*` | `basic.hh`, `basic_polynomial.hh` | 核心运算 |
| 单除数多项式除法 | `pair_vec_div`, `basic.hh:697` | 参考，但不能直接用于多除数约化 |
| 系数类型 ZZ、QQ、Zp | `number/` | QQ 为核心计算域 |
| 单项式整除判定 `is_divexact` | `monomial.hh:68` | Normal Form 的关键操作 |
| 单项式 GCD | `polynomial_gcd.hh:20` | 乘积准则判定 |
| 类型转换 ZZ↔QQ | `polynomial_convert.hh:63,118` | ZZ→QQ→ZZ 工作流 |
| 内容提取 `cont()` | `polynomial_gcd.hh:586` | ZZ 结果本原化 |
| 总度数 `degree()` | `polynomial.hh:40` | Sugar 度数初始化 |

### 1.2 缺失组件

| 组件 | 说明 |
|------|------|
| 单项式 LCM | 已有 `gcd`（取 min），缺对称的 `lcm`（取 max） |
| 多除数约化（Normal Form） | `pair_vec_div` 仅支持单除数，需新增多除数约化 |
| S-多项式 | 新增 |
| Buchberger 算法主体 | 临界对管理、选择策略、剪枝、互约 |

### 1.3 已有代码中的注意事项

- `basic_monomial::operator/` 不检查整除性，可产生负指数（使用前需 `is_divexact` 判定）
- QQ 无 `inv()` 方法，且无 `polynomial /= Tc` 运算符，首一化需逐系数除以 LC
- `pp()` 函数（polynomial_gcd.hh:694）仅限 `lex_` 序，通用序下本原化仍需手动 `cont()` + 除法
- `pair_vec_div` 是完整多项式长除法（单除数），不能直接复用于多除数约化

---

## 2. 外部实现调研

### 2.1 系统对比

| 系统 | 核心算法 | F4/F5 | Sugar | Gebauer-Möller | 系数域 | 开源 |
|------|---------|-------|-------|---------------|--------|------|
| **Singular** | Buchberger + Mora + slimgb | 否 | 是 | 是 | Z, Q, Z/pZ, 扩域 | GPL |
| **FLINT** | 朴素 Buchberger | 否 | 否 | 否 | Z, Z/nZ | LGPL |
| **Maple** | 编译 F4 + Buchberger | 是(F4) | — | — | Q, Z/pZ | 否 |
| **Mathematica** | Buchberger + Gröbner Walk | 否 | 未公开 | 未公开 | Q, 近似数, 有理函数 | 否 |
| **CoCoA** | Buchberger | 否 | 是 | 是 | Q, Z/pZ, twin-float | GPL |
| **SymPy** | 改进 Buchberger + F5B | 部分(F5B) | LCM 选择 | 是 | Z, Q, Z/pZ | BSD |
| **Magma** | F4 (Steel) | 是(F4) | — | — | 全 | 否 |
| **msolve** | F4 | 是(F4) | — | — | Q, Z/pZ | GPLv2 |
| **Groebner.jl** | F4 | 是(F4) | — | — | Z/pZ, Q | 开源 |

### 2.2 各系统要点

**Singular**（最广泛使用的开源 GB 引擎）：
- 核心为高度优化的 Buchberger + Sugar + Gebauer-Möller
- slimgb 变体：引入加权长度函数（项数 + ecart + 系数大小），保持中间多项式"苗条"
- Hilbert 驱动策略：先算 degrevlex GB 提取 Hilbert 函数，再引导目标序计算
- FGLM：零维理想的序转换算法
- `groebner` 命令自动启发式选择最优策略
- **不实现 F4/F5**，靠优化 Buchberger 达到高性能

**FLINT**：
- 仅朴素 Buchberger（代码中自述 "naive implementation"）
- 无 Sugar、无 Gebauer-Möller
- 提供资源限制版本（最大基大小、项数、系数位数）
- GB 支持极简，但多项式算术底层极优（packed exponents）

**Maple**：
- 默认使用编译 C 实现的 F4（限 Q 和 Z/pZ, p<2^31）
- F4 为 Monte-Carlo 概率算法（错误概率 <10^-18）
- 也保留 Buchberger 用于增量计算场景
- 集成 FGb 库

**Mathematica**：
- 优化 Buchberger + Gröbner Walk（lex 序时自动先算 degrevlex 再转换）
- 内部使用稀疏分布式多项式格式
- 内部细节未公开

**CoCoA**：
- Buchberger + Sugar + Gebauer-Möller 的经典实现
- 单项式（power-product）为独立一等对象 `PPMonoidElem`
- 最古老的 GB 专用系统之一（1987 年至今）

### 2.3 算法现状总结

**F4**（Faugère 1999）：当前实践中的高性能方案
- 核心思想：将多个 S-多项式约化转为稀疏矩阵行消元，批量处理
- 优势：缓存友好、可并行、摊薄开销
- 适用场景：大规模系统（>20 多项式、>5 变量），特别是有限域上
- 实现代价：需稀疏矩阵、全局单项式哈希表、符号预处理，工程量大

**签名算法（F5 及后继）**：理论前沿
- 用"签名"标记跟踪计算历史，在约化前检测冗余
- 实际可用实现较少（SymPy 的 F5B 是少数公开实现之一）

**Buchberger + Sugar + Gebauer-Möller**：经典可靠方案
- Singular、CoCoA、Mathematica、SymPy 均采用或以此为基础
- Sugar 是单一最有效的优化（Giovini et al. 1991："always the best choice or very near to it"）
- Gebauer-Möller 准则在实践中通常消除 80-95% 的冗余对
- 对通用系数域（非有限域/Q）最自然
- 内存友好：逐个处理 S-多项式

---

## 3. 方案评估

### 3.1 Buchberger + Sugar + Gebauer-Möller（推荐）

**优点：**
- 实现复杂度适中，CLPoly 已有基础设施完全覆盖所需前置条件
- 与 Singular（最成功的开源 GB 引擎）采用相同的核心方案
- 对 CLPoly 的模板化序设计天然适配
- 调试验证直接：约化 GB 唯一，可与 Mathematica 逐项对比
- 支持增量计算、任意系数域、任意序（含未来可能的局部序）

**已知局限：**
- 大规模系统（>20 多项式、>5 变量、高次）性能不及 F4
- QQ 上中间系数膨胀（可后续用模方法优化）

**未来演进路径（按优先级）：**
1. 模方法：Z/pZ 计算 + CRT 重构 Q 上结果（最大实用加速）
2. FGLM：零维理想 degrevlex → lex 序转换
3. Hilbert 驱动：利用 Hilbert 函数引导对选择
4. F4：稀疏矩阵 + 符号预处理（工程量大，性能质变）

### 3.2 F4（备选，不推荐作为首次实现）

**优点：** 大规模问题性能优越
**缺点：** 实现工程量大（稀疏矩阵、哈希表、符号预处理），与 CLPoly 现有数据结构差距较大

### 3.3 结论

**推荐 Buchberger + Sugar + Gebauer-Möller 作为首次实现。** 这是 Singular、CoCoA、Mathematica 等主流系统采用的经典方案，与 CLPoly 的模板化设计和已有基础设施高度契合。后续可按演进路径逐步优化。

---

## 4. 参考文献

- Buchberger, "Ein Algorithmus zum Auffinden der Basiselemente des Restklassenringes nach einem nulldimensionalen Polynomideal", PhD thesis, 1965
- Giovini, Mora, Niesi, Robbiano, Traverso, "One sugar cube, please", ISSAC 1991
- Gebauer & Möller, "On an Installation of Buchberger's Algorithm", J. Symbolic Comput. 1988
- Faugère, "A new efficient algorithm for computing Gröbner bases (F4)", J. Pure Applied Algebra 1999
- Eder & Faugère, "A survey on signature-based algorithms for computing Gröbner bases", J. Symbolic Comput. 2017
- Brickenstein, "Slimgb: Gröbner bases with slim polynomials", Rev. Mat. Complut. 2010
- Cox, Little & O'Shea, "Ideals, Varieties, and Algorithms", Springer
- Becker & Weispfenning, "Gröbner Bases: A Computational Approach to Commutative Algebra", Springer GTM 141, 1993
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra", Kluwer 1992
- von zur Gathen & Gerhard, "Modern Computer Algebra", Cambridge 2013
