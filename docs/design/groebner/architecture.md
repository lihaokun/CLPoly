# Gröbner 基架构文档

> 前置：[调研报告](research.md)
> 算法方案：Buchberger + Sugar 策略 + Gebauer-Möller 剪枝

---

## 1. 核心流程

### 1.1 总体数据流

```
用户输入: vector<polynomial> generators
        │
        ▼
   ┌─────────────┐
   │ 入口分派     │  groebner_basis(generators)
   │ ZZ → QQ 转换 │  ZZ 输入先转 QQ
   └──────┬──────┘
          │  vector<polynomial_QQ>
          ▼
   ┌─────────────┐
   │ Buchberger  │  __buchberger(F)
   │ 主循环       │
   └──────┬──────┘
          │  vector<polynomial_QQ>（约化 GB）
          ▼
   ┌─────────────┐
   │ 出口转换     │  QQ → ZZ 本原化（如需）
   └──────┬──────┘
          │
          ▼
用户输出: vector<polynomial>（约化 Gröbner 基）
```

### 1.2 Buchberger 主循环内部流程

```
输入 F（非零、首一化）
        │
        ▼
   初始化 sugar 度数
   初始化临界对集合 Pairs（Gebauer-Möller 过滤）
        │
        ▼
   ┌──────────────────────────────┐
   │ while Pairs 非空:            │
   │   1. 选对: 取 (sugar, lcm_deg) 最小的对 (i,j)
   │   2. S-多项式: h ← S(G[i], G[j])
   │   3. 约化: h ← NF(h, G)（带 sugar 跟踪）
   │   4. 若 h ≠ 0:
   │      a. 首一化 h
   │      b. Update: 用 h 更新 G 和 Pairs（Gebauer-Möller）
   └──────────────────────────────┘
        │
        ▼
   互约 __interreduce(G)
        │
        ▼
   输出: 约化 Gröbner 基
```

---

## 2. 模块划分

共 6 个模块，依赖关系如下（箭头表示"依赖于"）：

```
  M6: API 层 ──→ M5, M3
       │
       ▼
  M5: Buchberger 主循环 ──→ M2, M3, M4
       │
       ├──→ M2: S-多项式 ──→ M1
       ├──→ M3: Normal Form（无外部依赖）
       └──→ M4: 临界对管理 ──→ M1
                                │
                                ▼
                          M1: 单项式 LCM（底层，无依赖）
```

层次：
- 底层：M1（单项式 LCM）、M3（Normal Form）—— 无模块间依赖
- 中层：M2（依赖 M1）、M4（依赖 M1）
- 核心层：M5（依赖 M2、M3、M4）
- 顶层：M6（依赖 M5、M3）

---

## 3. 模块功能规约

### M1: 单项式 LCM

```
模块名称：单项式 LCM

功能描述：计算两个单项式的最小公倍式（逐变量取指数最大值）

前置条件（Requires）：
  - m1, m2 是同一 compare 类型的 basic_monomial

后置条件（Ensures）：
  - 返回 L = lcm(m1, m2)
  - m1 | L 且 m2 | L（即 L 的每个变量指数 ≥ m1 和 m2 中对应指数）
  - L 是满足上述条件的最小单项式（即 deg(L) = Σ max(e1_i, e2_i)）

副作用：无
```

**所在文件：** `clpoly/monomial.hh`（与已有 `gcd` 对称）

### M2: S-多项式

```
模块名称：S-多项式

功能描述：计算两个多项式的 S-多项式

前置条件（Requires）：
  - f, g 非零多项式
  - Tc 为域（QQ 或 Zp），保证系数精确除法

后置条件（Ensures）：
  - 返回 S(f, g) = (L/LT(f))·f - (L/LT(g))·g
  - 其中 L = lcm(LM(f), LM(g))
  - S(f, g) = 0，或 LM(S(f,g)) < L（首项被消去）

副作用：无
```

### M3: Normal Form（多除数约化）

```
模块名称：多除数约化

功能描述：计算多项式 f 关于多项式集合 G 的 Normal Form

前置条件（Requires）：
  - G 中各元素非零
  - Tc 为域，保证系数精确除法

后置条件（Ensures）：
  - 返回 r = NF(f, G)
  - f - r ∈ <G>（f 与 r 模理想 <G> 同余）
  - r 的任何单项式都不被 G 中任何元素的首项单项式整除
  - 特别地：若 f ∈ <G> 且 G 是 Gröbner 基，则 r = 0

副作用：无

变体：
  - 基本版 __normal_form(f, G)：返回 NF（供 M6 Ideal 类和 M5 互约使用）
  - Sugar 跟踪版 __normal_form_with_sugar(f, G, sugar, f_sugar)：
    同时更新 f 的 sugar 度数（供 M5 主循环使用）
```

### M4: 临界对管理

```
模块名称：临界对管理（Sugar 策略 + Gebauer-Möller 剪枝）

功能描述：维护临界对集合，提供选择和更新操作

前置条件（Requires）：
  - G 为当前基，sugar 为各元素的 sugar 度数
  - |sugar| = |G|

后置条件（Ensures）：
  - 选择操作：返回并移除 (sugar, lcm_deg) 最小的临界对
  - 更新操作（新元素 h 加入时）：
    a. 用 Gebauer-Möller 三准则过滤，去除可证明约化到零的对：
       - 乘积准则：gcd(LM(f), LM(g)) = 1 则 S(f,g) →_G 0
       - 链准则：存在 k 使 LM(G[k]) | lcm(LM(G[i]), LM(G[j]))
         且 (i,k) 和 (k,j) 已被或将被处理，则 (i,j) 可丢弃
       - LCM 准则：LM(h) | lcm(LM(G[i]), LM(G[j])) 且两个子 lcm
         都严格更小，则旧对 (i,j) 可丢弃
    b. 新对的 sugar 度数正确计算
    c. 旧对中被 h 淘汰的对被移除

不变式（Invariants）：
  - 对集合中不包含可被 Gebauer-Möller 准则淘汰的对
  - 每个对的 sugar 和 lcm_deg 字段与 G 中对应元素一致

副作用：修改 G（追加 h）、修改 Pairs、修改 sugar（追加 h_sugar）
```

### M5: Buchberger 主循环

```
模块名称：Buchberger 主循环

功能描述：实现 Buchberger 算法，计算约化 Gröbner 基

前置条件（Requires）：
  - F 中各元素为 QQ（或域）上的多项式
  - F 可包含零多项式（由 M5 内部过滤，调用方无需预处理）

后置条件（Ensures）：
  - 返回 G = 约化 Gröbner 基
  - <G> = <F>（生成相同理想）
  - G 中每个元素首一（LC = 1）
  - G 中无冗余生成元（无 LM(gi) | LM(gj), i ≠ j）
  - G 中每个元素关于 G\{g} 完全约化
  - G 是唯一的（约化 GB 对给定理想和序唯一）

副作用：无（F 按值传入）

子流程：
  - 互约 __interreduce(G)：将 GB 化简为约化 GB
```

### M6: API 层

```
模块名称：API 层

功能描述：提供用户入口，处理类型转换，封装 Ideal 类

前置条件（Requires）：
  - 输入为 QQ 或 ZZ 上的多项式集合

后置条件（Ensures）：
  - QQ 入口：直接返回 QQ 上的约化 GB
  - ZZ 入口：返回 ZZ 上的约化 GB，各元素本原且 lc > 0
  - Ideal 类：惰性计算 GB，缓存结果
    - contains(f) 返回 true ⟺ f ∈ <generators>
    - reduce(f) 返回 NF(f, GB)

副作用：Ideal 类内部修改缓存（mutable）
```

---

## 4. 模块间接口规约

### 接口 I1: M1 → M2（LCM 供 S-多项式使用）

```
接口：M1(单项式 LCM) → M2(S-多项式)

输入数据：两个 basic_monomial<comp>（即 LM(f) 和 LM(g)）
输出数据：basic_monomial<comp>（即 lcm(LM(f), LM(g))）

协议约定：
  - 调用方：保证 m1, m2 为合法的非零单项式
  - 被调用方：保证返回值的度数正确（deg = Σ max(e1_i, e2_i)），无需 re_deg()
```

### 接口 I2: M1 → M4（LCM 供临界对管理使用）

```
接口：M1(单项式 LCM) → M4(临界对管理)

输入数据：两个 basic_monomial<comp>
输出数据：basic_monomial<comp>

协议约定：
  - 调用方：用于计算 lcm_deg 和 Gebauer-Möller 准则判定
  - 被调用方：同 I1
```

### 接口 I3: M2 → M5（S-多项式供主循环使用）

```
接口：M2(S-多项式) → M5(Buchberger 主循环)

输入数据：G 中两个多项式 G[i], G[j]
输出数据：polynomial_<Tc, comp>（S-多项式）

协议约定：
  - 调用方：保证 G[i], G[j] 非零
  - 被调用方：返回 S(G[i], G[j])，首项已消去
```

### 接口 I4: M3 → M5（Normal Form sugar 版供主循环使用）

```
接口：M3(__normal_form_with_sugar) → M5(Buchberger 主循环)

输入数据：
  - h: 待约化多项式（S-多项式）
  - G: 当前基
  - sugar: 各基元素的 sugar 度数
  - h_sugar: h 的 sugar 度数（引用，可被修改）
输出数据：polynomial_<Tc, comp>（约化结果）

协议约定：
  - 调用方：保证 G 中各元素非零，|sugar| = |G|
  - 被调用方：
    a. 返回 r = NF(h, G)
    b. h_sugar 更新为约化过程中的最终 sugar 度数
    c. 约化规则：sugar(h) ← max(sugar(h), sugar(g) + deg(LM(h)) - deg(LM(g)))
```

### 接口 I5a: M4 → M5（选择操作）

```
接口：M4(临界对管理) → M5(Buchberger 主循环) —— 选择

输入数据：当前对集合 Pairs（有序集合，细化阶段选用 vector + 线性扫描实现）
输出数据：__critical_pair（(sugar, lcm_deg) 最小的对，从集合中移除）

协议约定：
  - 调用方：保证 Pairs 非空
  - 被调用方：返回并移除优先级最高（sugar 最小，sugar 相同则 lcm_deg 最小）的对
```

### 接口 I5b: M4 → M5（更新操作）

```
接口：M4(临界对管理) → M5(Buchberger 主循环) —— 更新

输入数据：G, Pairs, h（新元素）, sugar, h_sugar
输出数据：修改后的 G（追加 h）, Pairs（增删对）, sugar（追加 h_sugar）

协议约定：
  - 调用方：保证 h 非零且已首一化
  - 被调用方：
    a. 用 Gebauer-Möller 三准则过滤
    b. 将 h 追加到 G 末尾
    c. 新对的 sugar 正确计算
```

### 接口 I6: M5 → M6（Buchberger 供 API 层使用）

```
接口：M5(Buchberger) → M6(API 层)

输入数据：vector<polynomial_<QQ, comp>>（生成元，QQ 系数）
输出数据：vector<polynomial_<QQ, comp>>（约化 GB，首一化）

协议约定：
  - 调用方（API 层）：
    a. QQ 入口：直接传入
    b. ZZ 入口：先用 poly_convert 转为 QQ
  - 被调用方：返回约化 GB，满足 M5 功能规约中的所有后置条件
```

### 接口 I7: M3 → M6（Normal Form 基本版供 Ideal 类使用）

```
接口：M3(__normal_form) → M6(Ideal 类)

输入数据：f（待约化多项式）, G（约化 GB，由 Ideal 缓存）
输出数据：polynomial_<Tc, comp>（NF）

协议约定：
  - 调用方：保证 G 是约化 GB（通过 Ideal 内部缓存保证）
  - 被调用方：返回 NF(f, G)，满足 M3 功能规约（基本版，不跟踪 sugar）
```

---

## 5. 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 核心算法 | Buchberger | 经典方案，Singular/CoCoA/Mathematica 均采用 |
| 选择策略 | Sugar + 度数 | "always the best choice"（Giovini et al. 1991） |
| 剪枝 | Gebauer-Möller 三准则 | 实践中消除 80-95% 冗余对 |
| 计算域 | QQ 上实现核心 | 域上运算精确，NF 无需伪除法 |
| ZZ 支持 | ZZ→QQ→ZZ | 复用已有 poly_convert |
| 单项式序 | 模板化 comp | 同一实现支持 lex、grlex 及任意序 |
| 文件组织 | 新建 `groebner.hh` | 与 polynomial_gcd.hh、polynomial_factorize.hh 并列 |
| 内部函数 | `__` 前缀 | 与 `__polynomial_GCD` 风格一致 |
| 首一化 | `__make_monic(f)` | 无 `polynomial /= Tc`，用辅助函数逐系数除以 LC |
| LCM 位置 | `monomial.hh` | 与已有 gcd(basic_monomial) 对称 |

---

## 6. 文件组织

### 6.1 新增文件

| 文件 | 内容 |
|------|------|
| `clpoly/groebner.hh` | M2-M6 全部代码 |
| `test/test_groebner.cc` | 测试 |

### 6.2 修改文件

| 文件 | 修改 |
|------|------|
| `clpoly/monomial.hh` | 新增 M1（lcm） |
| `clpoly/clpoly.hh` | 新增 `#include "groebner.hh"` |
| `test/run_all_tests.sh` | 注册 test_groebner |
| `.gitignore` | 注册 test/test_groebner |
| `README.md` | Gröbner 基 → 支持 |

### 6.3 groebner.hh 内部组织

```
头部（pragma, include, namespace）
  │
  ├── M2: __s_polynomial
  ├── M3: __normal_form, __normal_form_with_sugar
  ├── M4: __critical_pair, __update
  ├── M5: __interreduce, __buchberger
  └── M6: s_polynomial, normal_form, groebner_basis, Ideal
```
