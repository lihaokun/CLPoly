# 多项式系统对比研究与 CLPoly 发展路线

> 本文档归档了对 Maple、Mathematica、FLINT 及 CLPoly 四个系统在多项式表示与运算上的深入对比研究，
> 并据此制定 CLPoly 的短期与长期发展计划。

---

## 目录

1. [设计哲学总览](#1-设计哲学总览)
2. [数据结构对比](#2-数据结构对比)
3. [单项式编码](#3-单项式编码)
4. [核心算法对比](#4-核心算法对比)
5. [各系统独特机制](#5-各系统独特机制)
6. [性能基准参考](#6-性能基准参考)
7. [CLPoly 现状分析](#7-clpoly-现状分析)
8. [短期计划：因式分解与测试](#8-短期计划因式分解与测试)
9. [长期计划：新一代架构](#9-长期计划新一代架构)
10. [关键参考文献](#10-关键参考文献)

---

## 1. 设计哲学总览

| | Maple | Mathematica | CLPoly | FLINT |
|--|-------|-------------|--------|-------|
| 语言 | C 内核 + Maple 库 | C 内核 + WL 库 | C++ | C |
| 多项式是否专用类型 | **是**（POLY，v17+） | **否**（通用表达式树） | **是** | **是** |
| 设计取向 | 务实——为性能引入专用类型 | 统一——一切皆表达式 | 泛型——C++ 模板 | 极致性能 |
| 内存布局 | AoS（交织） | 指针树 | AoS（交织） | SoA（分离） |

**核心观察**：Maple 和 FLINT 代表了性能优先路线，Mathematica 代表了灵活性优先路线。
CLPoly 的长期目标是在两者之间找到最优平衡点。

---

## 2. 数据结构对比

### 2.1 Maple：三代表示的演进

#### 第一代（v1-13）—— Sum-of-Products DAG

所有 Maple 对象存储为 DAG（有向无环图），多项式用 SUM 和 PROD 两种节点：

```
SUM 节点:  [header | mono_1 | coeff_1 | mono_2 | coeff_2 | ... | const]
PROD 节点: [header | base_1 | exp_1  | base_2 | exp_2  | ...]

例: 3x²y + 5z - 7
SUM [ptr→(PROD x 2 y 1)] [3] [ptr→(PROD z 1)] [5] [1] [-7]
```

关键设计：
- 全局**唯一性哈希表（simplification table）**：相同对象只存一份
- 等价性测试 = 指针比较 = **O(1)**
- 交换哈希函数：`hash(2x + 3y) == hash(3y + 2x)`

缺点：
- 单项式运算需要分配内存 + 查找哈希表
- **占总计算时间 >70%**

#### 第二代（v14-16）—— SDMP 库（转换模式）

引入位压缩 SDMP（Sparse Distributed Multivariate Polynomial）格式：

```
用户输入 → Sum-of-Products → 转 SDMP → 快速运算 → 转回 Sum-of-Products
```

- 除法速度提升 **250x**
- 但 SoP ↔ SDMP 转换开销仍是瓶颈

#### 第三代（v17+）—— POLY 内核类型

```
[header] [ptr→变量列表] [packed_mono_1 coeff_1] [packed_mono_2 coeff_2] ...
```

- 位压缩单项式：3 变量时每个域 16 bits，单项式比较 = 1 条机器指令
- 内联小整数系数：`|x| < 2^62` → 存 `2x+1`（奇数标记），与指针（偶数）区分
- 无需格式转换，比 SoP 紧凑 **4x**
- 可处理**数亿项**的多项式

### 2.2 Mathematica —— 通用表达式树

```mathematica
3x²y + 5z - 7
↓ FullForm
Plus[Times[3, Power[x, 2], y], Times[5, z], -7]
```

内部结构：
```
Plus ─┬─ Times ─┬─ 3
      │         ├─ Power ─┬─ x
      │         │         └─ 2
      │         └─ y
      ├─ Times ─┬─ 5
      │         └─ z
      └─ -7
```

- 每个节点 = 连续指针数组 + 哈希码 + 引用计数
- `Plus` 和 `Times` 参数**自动排序**为规范序（`Orderless` 属性）
- 没有专用多项式类型——`PolynomialQ` 只是语法检查
- 有未公开的 `Internal'DistributedTermsList` 但非主要路径

### 2.3 CLPoly —— 泛型 C++ 模板

```cpp
// basic_polynomial<compare_type, T> = vector<pair<compare_type, T>>
// 多项式 = vector<pair<monomial, coefficient>>
// monomial = vector<pair<variable, degree>>
```

- 动态变量集合：不同多项式可有不同变量
- `pair_vec` 有序对向量作为统一容器
- 乘法时可临时压缩为位压缩格式

### 2.4 FLINT —— 极致性能 C 实现

```c
// 三级表示：mpoly（通用）→ 内部选择
// 存储：系数数组 + 指数数组（SoA 分离）
// 指数：位压缩到 ulong 数组
```

- Structure of Arrays（SoA）布局，缓存友好
- 自动 repack：溢出时扩展位宽而非回退

---

## 3. 单项式编码

### 同一单项式 x²yz³ 在四个系统中的表示

```
Maple POLY:   0x0006000200010003   ← 一个 uint64
              [deg=6|x=2|y=1|z=3]    每域 16 bits

FLINT:        0x0006000200010003   ← 一个或多个 ulong
              [deg=6|x=2|y=1|z=3]    可 repack 到更宽

CLPoly:       vector<pair> = [(x,2),(y,1),(z,3)]
              __deg = 6              变长符号化

Mathematica:  Times[Power[x,2], y, Power[z,3]]
              ← 表达式树，多个堆对象
```

### 操作开销对比

| 操作 | Maple POLY | FLINT | CLPoly | Mathematica |
|------|-----------|-------|--------|-------------|
| 比较 | **1 条指令** | **1 条指令** | O(k) 遍历 | O(k) 树遍历 + 哈希 |
| 乘法（指数加） | **1 条指令** | **1 条指令** | O(k) 逐对加 | 构造新 Times 节点 |
| 内存/项 | ~16 bytes | ~16 bytes | ~80+ bytes | ~200+ bytes |
| 溢出处理 | 回退 SoP | repack 到更宽 | 无 | 不适用 |

### 位压缩容量（64 位机器）

| 变量数 (n) | 每指数位数 | 最大次数 |
|-----------|-----------|---------|
| 1 | 32 | 4,294,967,295 |
| 2 | 21 | 2,097,151 |
| 3 | 16 | 65,535 |
| 5 | 10 | 1,023 |
| 7 | 8 | 255 |
| 15 | 4 | 15 |

---

## 4. 核心算法对比

### 4.1 加法

四个系统都使用 **O(n+m) 双指针归并**，但每步开销差异巨大：

```
一次"比较 + 前进"的开销：

Maple POLY:   cmp(uint64, uint64)              ≈ 1 ns
FLINT:        cmp(ulong, ulong)                ≈ 1 ns
CLPoly:       遍历 vector<pair> 逐对比较        ≈ 20-50 ns
Mathematica:  表达式树遍历 + 哈希比较 + 求值器  ≈ 100-500 ns
```

### 4.2 乘法

#### 算法选择策略

```
Maple (v18+):
  ┌─ 密集? → Kronecker 替换 → GMP 大整数乘法 (FFT/Karatsuba)
  ├─ 短底数高次幂? → 专用幂展开算法
  └─ 稀疏? → Monagan-Pearce 堆算法 (多线程)

FLINT:
  ┌─ 尝试 Array 乘法 → 成功则返回
  ├─ 尝试 Dense 乘法 (Kronecker) → 成功则返回
  └─ Heap 乘法 (Monagan-Pearce) → 始终成功 (多线程)

CLPoly:
  ┌─ 能位压缩? → 压缩后堆乘法
  └─ 不能? → 普通堆乘法

Mathematica:
  ┌─ 系数列表? → ListConvolve (FFT / 二进制分段)
  └─ 符号表达式? → Expand (表达式树操作)
```

#### 堆乘法实现关键差异

| | 经典 Johnson | Monagan-Pearce 改进 |
|--|-------------|-------------------|
| 比较次数 | O(nm·log min(n,m)) | **O(nm)** |
| 堆存储 | 完整乘积项 | **仅存索引 (i,j)** |
| 内存分配 | 频繁 | **预分配，零垃圾** |
| 使用者 | CLPoly | **Maple, FLINT** |

CLPoly 当前使用传统堆算法，堆中存完整 `VHC` 结构体；Maple/FLINT 仅存索引对，实际乘积在提取时才计算。

#### 密集乘法路径

| | 方法 | 大整数库 |
|--|------|---------|
| Maple | Kronecker 替换 | GMP |
| FLINT | Kronecker 替换 / Array 累加 | GMP |
| CLPoly | 无密集路径 | — |
| Mathematica | `ListConvolve` + 二进制分段 | GMP (v5.0+) |

### 4.3 除法

| | 算法 | 能力 |
|--|------|------|
| **Maple** | Monagan-Pearce 堆除法 | 完整多项式除法，O(#q·#g·log #q) |
| **FLINT** | 堆除法 + Array 回退 | 完整多项式除法 + 带余除法 |
| **CLPoly** | `pair_vec_div` | 仅**单项式除法**（除数只有一项） |
| **Mathematica** | 内核实现 | `PolynomialQuotientRemainder` |

### 4.4 堆 vs Geobuckets

| 属性 | 堆 (Maple/FLINT) | Geobuckets (Singular) |
|------|-----------------|----------------------|
| 单项式比较数 | 较多 | **较少**（约一半） |
| 内存使用 | **少得多** | 较多 |
| 缓存性能 | **更好**（L1/L2） | 需 L2 容纳 |
| 大规模稀疏 | **更快** | 较慢（缓存不命中） |
| 密集情形 | 相当 | 相当 |

---

## 5. 各系统独特机制

### 5.1 Maple：唯一性哈希表

```
任何对象创建后：
  1. 递归简化
  2. 计算交换哈希（hash(a+b) = hash(b+a)）
  3. 查全局哈希表 → 如已存在，返回现有指针

结果：相等性测试 = 指针比较 = O(1)
```

POLY 表示下仍保留——整个 POLY 对象也进入简化表。

### 5.2 Mathematica：表达式求值器

```
每个表达式 h[e1,e2,...] 的求值：
  1. 求值 head h
  2. 求值每个 ei
  3. 应用 Flat、Orderless 等属性
  4. 应用用户定义的 DownValues/UpValues
  5. 应用内建规则
  6. 重复直到不变
```

多项式每步运算结果都经完整求值器——这是 Mathematica 在纯多项式运算上慢的根本原因，
但也赋予了它无与伦比的符号规则匹配能力。

### 5.3 CLPoly：动态变量系统

```cpp
// 不同多项式可以有不同变量集合
poly1: {x, y}    → 3x²y + 1
poly2: {y, z}    → 2yz - 5
poly1 * poly2:   → 自动合并为 {x, y, z}
```

其他三个系统都需要预先定义变量集合或统一变量环境。
这是 CLPoly 在灵活性上的独特优势。

### 5.4 FLINT：溢出位 + repack

```
[0|1100100]  +  [0|0011110]  =  [1|0000010]
 正常(100)      正常(30)        溢出！→ repack 到更宽
```

Maple 溢出时回退到 Sum-of-Products；FLINT 则在位压缩框架内解决。

---

## 6. 性能基准参考

### Fateman 基准

`f1 = (1+x+y+z)^20 + 1`，计算 `f1 × (f1+1)`

| 系统 | 大致相对速度 | 说明 |
|------|-------------|------|
| **FLINT** | 1x（最快） | 三级算法选择 + 位压缩 |
| **Maple 17+** | ~1-2x | POLY + Monagan-Pearce 堆 |
| **CLPoly** | ~5-20x（估计） | 临时压缩 + 传统堆 |
| **Mathematica** | ~75x（Monagan 报告） | 表达式树开销 |

### 性能-灵活性轴

```
性能轴（快→慢）:  FLINT ≈ Maple POLY >> CLPoly >> Mathematica
灵活轴（高→低）:  Mathematica > CLPoly > Maple > FLINT
```

---

## 7. CLPoly 现状分析

### 7.1 已有能力

| 模块 | 文件 | 状态 |
|------|------|------|
| 多项式基本运算 | `polynomial.hh`, `basic_polynomial.hh` | 完善 |
| 单项式运算与序 | `monomial.hh`, `monomial_order.hh` | 完善 |
| GCD 计算 | `polynomial_gcd.hh` (31K，最大模块) | 完善 |
| 子结式 / 结式 | `resultant.hh` | 完善 |
| 实根隔离 | `realroot.hh` | 完善 |
| 区间运算 | `interval.hh` | 完善 |
| 特征列 | `charset.hh` | 完善 |
| 一元多项式 | `upolynomial.hh` | 基本 |
| 随机多项式生成 | `random.hh` | 基本 |

### 7.2 缺失能力

> ⚠️ **[部分过时 2026-02-23]** 因式分解和 Hensel 提升已实现，见下表更新。

| 能力 | 状态 | 对标系统 |
|------|------|---------|
| **因式分解** | ✅ **已实现**（单变量 van Hoeij LLL + 多变量 Wang EEZ） | Maple `factor()`, Mathematica `Factor[]` |
| **完整多项式除法** | 仅单项式除法 | Maple/FLINT 的堆除法 |
| **密集乘法路径** | 缺失 | Kronecker 替换 |
| **位压缩默认存储** | 仅临时使用 | Maple POLY, FLINT |
| **并行运算** | 缺失 | Maple 多线程堆乘法 |
| **Hensel 提升** | ✅ **已实现**（二次提升 + 两阶段精度方案） | 因式分解核心组件 |

### 7.3 与 Maple 演进的类比

CLPoly 当前状态类似 **Maple v14-16**：
- 知道位压缩是正确方向（已有临时压缩）
- 但只在局部应用，未提升为默认存储格式
- 核心算法可用但未达到最优

---

## 8. 短期计划：因式分解与测试

### 8.1 因式分解实现路径

多元多项式因式分解是分层构建的，需要按以下顺序实现：

```
Layer 1: 一元整数多项式因式分解 (Z[x])
  ├─ 无平方分解 (square-free factorization) ← CLPoly 已有基础
  ├─ Berlekamp / Cantor-Zassenhaus 算法 (Fp[x])
  ├─ Hensel 提升 (Fp[x] → Z[x])
  └─ Zassenhaus 组合 / van Hoeij LLL 格基约化

Layer 2: 一元有理多项式因式分解 (Q[x])
  ├─ 去分母 → 归结到 Z[x]
  └─ 内容提取 + 本原化

Layer 3: 多元多项式因式分解 (Z[x1,...,xn])
  ├─ 赋值齐次化 → 归结到低变量数
  ├─ Wang's EEZ (Extended Ezgcd) 算法
  │   或 Musser-Kaltofen-Trager 多元提升
  └─ Leading coefficient correction

Layer 4: 有理系数多元因式分解 (Q[x1,...,xn])
  └─ 去分母 + Layer 3
```

### 8.2 具体实施步骤

#### 阶段一：基础设施（~2 周）

1. **完善 `upolynomial` 模块**
   - 实现完整的一元多项式除法 (`divmod`)
   - 实现模素数运算 (`Fp[x]` 基本运算)
   - 实现一元 GCD over Fp

2. **实现无平方分解**
   - `squarefree_factorization(f)` → 返回 `[(f1, 1), (f2, 2), ...]`
   - CLPoly 已有 `polynomial_squarefree_test.cc`，需扩展为完整算法

#### 阶段二：Fp[x] 因式分解（~2 周）

3. **Cantor-Zassenhaus 算法**
   - 等次分解 (distinct-degree factorization)
   - 等次分裂 (equal-degree splitting)
   - 对奇特征 p 和 p=2 分别处理

#### 阶段三：Z[x] 因式分解（~3 周）

4. **Hensel 提升**
   - 线性提升 (linear Hensel lifting)
   - 二次提升 (quadratic Hensel lifting, 更快)

5. **因子组合**
   - Zassenhaus 枚举法（小情形）
   - van Hoeij LLL 算法（通用情形，避免指数爆炸）

#### 阶段四：多元提升（~4 周）

6. **多元因式分解**
   - 赋值策略：选好的赋值点使一元因式分解无粘连
   - 多元 Hensel 提升 (Wang's algorithm)
   - 首项系数修正 (leading coefficient correction)

### 8.3 测试计划

#### 现有测试覆盖

```
test/
├── polynomial_squarefree_test.cc   (722 bytes, 待扩充)
├── polynomial_GCD_test.cc          (3.5K)
├── polynomial_test.cc              (8.7K)
├── upolynomial_test.cc             (4.3K)
└── ... (共 20 个测试文件)
```

#### 需要新增的测试

| 测试文件 | 覆盖内容 |
|---------|---------|
| `upolynomial_divmod_test.cc` | 一元除法的正确性和边界情况 |
| `upolynomial_modp_test.cc` | Fp[x] 运算，不同素数 p |
| `factorization_univariate_test.cc` | Z[x] 因式分解正确性 |
| `factorization_multivariate_test.cc` | Z[x1,...,xn] 因式分解 |
| `hensel_lifting_test.cc` | Hensel 提升过程验证 |
| `factorization_benchmark.cc` | 性能基准（对标 Maple/Mathematica 结果） |

#### 测试用例来源

- Mathematica notebook `test/未命名-1.nb` 中已有的测试数据
- Maple 文档中的经典因式分解示例
- FLINT 测试套件中的因式分解用例
- Swinnerton-Dyer 多项式（因式分解的经典压力测试）
- Lenstra 的 LLL 论文中的示例

---

## 9. 长期计划：新一代架构

### 9.1 目标

> 兼顾 Mathematica 的灵活性与 Maple 的速度

这意味着：
- **用户面**：像 Mathematica 一样自然地处理符号表达式
- **内核面**：像 Maple POLY 一样高效地存储和运算多项式

### 9.2 双层架构设想

```
┌────────────────────────────────────────────────┐
│           用户层 (Expression Layer)              │
│                                                  │
│  • 通用表达式树（类似 Mathematica）                │
│  • 支持 sin, exp, 矩阵, 微分方程 等               │
│  • 动态变量系统（保留 CLPoly 现有优势）             │
│  • 模式匹配和规则替换                              │
│                                                  │
├────────────────────────────────────────────────┤
│         多项式加速层 (Polynomial Fast Path)       │
│                                                  │
│  • 检测到纯多项式时自动下沉到快速路径               │
│  • 位压缩存储（Maple POLY 风格）                   │
│  • Monagan-Pearce 堆乘法                          │
│  • 堆除法                                        │
│  • Kronecker 替换（密集路径）                      │
│  • 结果自动提升回表达式层                          │
│                                                  │
├────────────────────────────────────────────────┤
│           底层 (Arithmetic Backend)              │
│                                                  │
│  • GMP 大整数                                     │
│  • 机器字内联小整数（Maple 的奇偶标记法）           │
│  • SIMD 加速（AVX2/AVX-512 向量化单项式运算）      │
│  • 多线程支持                                     │
│                                                  │
└────────────────────────────────────────────────┘
```

### 9.3 潜在技术路径

#### 路径 A：渐进式改造（推荐）

从现有 CLPoly 代码库出发，逐步引入高性能组件：

```
Phase 1: 将位压缩提升为默认存储
  └─ 这正是 Maple v14 → v17 所走的路
  └─ 保留 vector<pair> 作为回退路径
  └─ 预期收益：运算速度 5-10x 提升

Phase 2: 引入 Monagan-Pearce 堆算法
  └─ 替换当前的传统堆实现
  └─ 堆中仅存索引对，预分配内存
  └─ 预期收益：乘法速度 2-3x 提升

Phase 3: 密集路径 (Kronecker 替换)
  └─ 自动检测密集多项式
  └─ 通过 GMP 的 FFT 乘法加速
  └─ 预期收益：密集情形 10-100x 提升

Phase 4: 表达式层
  └─ 在多项式层之上构建通用表达式
  └─ 多项式运算自动路由到快速路径

Phase 5: 并行化
  └─ Monagan-Pearce 并行堆乘法
  └─ 多核 GCD、因式分解
```

#### 路径 B：重新设计（高风险高回报）

从零设计新的内核：

```
• 统一值表示：tagged union = {SmallInt | BigInt | Polynomial | Expression | ...}
• 全局唯一性表（Maple 风格）+ 引用计数/GC
• JIT 编译热路径
• 风险：工程量巨大，可能 2-3 年才能达到现有功能覆盖
```

#### 路径 C：FLINT 后端 + CLPoly 前端

```
• CLPoly 保留现有 API 和灵活性
• 底层多项式运算委托给 FLINT
• 需要处理 CLPoly 动态变量 ↔ FLINT 固定变量环 的转换
• 风险：两个系统的数据模型差异可能导致大量转换开销
```

### 9.4 各路径评估

| | 路径 A（渐进） | 路径 B（重写） | 路径 C（FLINT 后端） |
|--|--------------|--------------|-------------------|
| 开发时间 | 中（6-12 月） | 长（2-3 年） | 短（3-6 月） |
| 风险 | **低** | 高 | 中 |
| 最终性能 | 高（~Maple 水平） | 极高（理论最优） | 极高（=FLINT） |
| 灵活性保持 | **完全保持** | 取决于设计 | 部分受限 |
| 代码可控性 | **完全自主** | 完全自主 | 依赖外部库 |

**建议**：采用路径 A 作为主路线，在 Phase 3 之后评估是否需要切换到路径 B。

### 9.5 里程碑规划

```
M1 (短期, 1-3 个月):
  ✓ 因式分解基本实现
  ✓ 测试套件完善
  ✓ 完整多项式除法

M2 (中期, 3-6 个月):
  □ 位压缩默认存储 (Phase 1)
  □ Monagan-Pearce 堆算法 (Phase 2)
  □ 性能基准套件

M3 (中期, 6-12 个月):
  □ Kronecker 替换 (Phase 3)
  □ 因式分解优化
  □ 性能达到 Maple 同一量级

M4 (长期, 12-24 个月):
  □ 表达式层设计与实现 (Phase 4)
  □ 并行化 (Phase 5)
  □ 符号积分 / 微分方程 等扩展
```

---

## 10. 关键参考文献

### Maple 多项式实现（Monagan & Pearce 系列）

1. **"Polynomial Division Using Dynamic Arrays, Heaps, and Packed Exponent Vectors"**
   CASC 2007, Springer LNCS 4770, pp. 295-315
   — 引入 SDMP 数据结构和位压缩方案

2. **"Parallel Sparse Polynomial Multiplication Using Heaps"**
   ISSAC 2009, ACM Press, pp. 263-270
   — 多核并行堆乘法，超线性加速

3. **"Sparse Polynomial Multiplication and Division in Maple 14"**
   ACM Communications in Computer Algebra, 44:4, 2010
   — 堆算法集成到 Maple expand/divide

4. **"Sparse Polynomial Division Using a Heap"**
   Journal of Symbolic Computation, 46:7, 2011, pp. 807-822
   — 无分数堆除法，与乘法同等复杂度

5. **"POLY: A New Polynomial Data Structure for Maple 17"**
   ACM Communications in Computer Algebra, 46:3/4, 2012
   — POLY 内核类型设计

6. **"Sparse Polynomial Powering Using Heaps"**
   CASC 2012, Springer LNCS 7442
   — 稀疏多项式幂运算

7. **"The Design of Maple's Sum-of-Products and POLY Data Structures"**
   ACM Communications in Computer Algebra, 48:3/4, 2014
   — 两代数据结构的权威描述

8. **"Speeding up Polynomial GCD, a Crucial Operation in Maple"**
   Maple Transactions, Vol. 2, No. 1, 2022
   — 模 GCD 算法 + POLY 系数编码细节

### 因式分解算法

9. **Zassenhaus, H.** "On Hensel Factorization I"
   Journal of Number Theory, 1969
   — 经典 Hensel 提升方法

10. **van Hoeij, M.** "Factoring Polynomials and the Knapsack Problem"
    Journal of Number Theory, 95:2, 2002
    — LLL 格基约化避免因子组合的指数爆炸

11. **Wang, P.S.** "An Improved Multivariate Polynomial Factoring Algorithm"
    Mathematics of Computation, 32:144, 1978
    — 多元 Hensel 提升 (EEZ 算法)

### 通用参考

12. **von zur Gathen, J. & Gerhard, J.** *Modern Computer Algebra*, 3rd ed., Cambridge, 2013
    — 计算机代数权威教材

13. **Geddes, K., Czapor, S., & Labahn, G.** *Algorithms for Computer Algebra*, Springer, 1992
    — 经典参考，涵盖因式分解完整理论

---

*文档创建日期: 2026-02-11*
*项目: CLPoly — C++ 多项式库*
*作者: CLPoly 开发团队*
