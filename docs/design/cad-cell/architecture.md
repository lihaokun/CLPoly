# CAD Cell 架构文档

## 1. 核心流程

```
输入: 多项式集 F ⊆ Z[x_1,...,x_n], 变量序 vars = [x_1,...,x_n]
                                      (x_1 lex最小 = 投影最先消去)

Phase 1: 投影（已实现）
  __project_full(F, vars) → allprojs[0..n-1]
    allprojs[k] = P_{k+1} ⊆ Z[x_1,...,x_{k+1}]
    所有多项式都含 vars[k]

Phase 2: 提升（本次实现）
  cad_tree 逐层构建:
    Level 0: 对 allprojs[n-1]（纯 x_n 的一元多项式）做 realroot
             → 根把 R 分为 Section/Sector 节点
    Level k (k≥1): 对 Level k-1 的每个叶节点:
             取采样点 α，代入 allprojs[n-1-k] 中的多项式
             → 得到关于 x_{n-k} 的一元多项式集
             → realroot → Section/Sector 子节点
    叶节点 = R^n 中的完整胞腔

输出: cad_tree 对象，可遍历所有叶节点（完整 cell）及其采样点
```

**变量层级映射**：

投影阶段 `allprojs` 的索引约定（`__project_full` 已实现）：
- `allprojs[0]` 含 `vars[0]`（lex 最小变量），多项式在 Z[x_1,...,x_n] 中
- `allprojs[k]` 含 `vars[k]`，多项式在 Z[x_1,...,x_{k+1},...] 中
- `allprojs[n-1]` 含 `vars[n-1]`（lex 最大变量），为一元多项式

提升阶段从 `allprojs[n-1]`（一元）开始向 `allprojs[0]` 推进：
- 树的 Level 0 对应 `allprojs[n-1]`（一元，无需代入）
- 树的 Level k 对应 `allprojs[n-1-k]`（需代入前 k 层的采样点）

## 2. 模块划分

### 2.1 `cad_tree` — 管理类（唯一对外接口）

```
类名: cad_tree<var_order>

模板参数:
  var_order — 变量比较器（同 lex_<var_order> 中的 var_order）

内部存储:
  _level_vars     : vector<variable>          — 提升序变量，_level_vars[k] = level k 对应的变量
  _level_polys    : vector<vector<poly_type>>  — 提升序投影多项式，_level_polys[k] = level k 的多项式
  _levels         : vector<vector<cad_node>>   — 分层节点，_levels[k] = level k 所有节点

  三个数组统一用提升序索引:
    level 0 = 一元多项式层（最先提升）
    level n-1 = 最后提升（叶层）
    _level_vars[k]、_level_polys[k]、_levels[k] 三者对应同一层

  构造时从外部输入转换:
    输入 vars[i] 按字典序（vars[0] = lex 最小）
    输入 allprojs[i] 与 vars 同序
    _level_vars[k] = vars[n-1-k]    （反转）
    _level_polys[k] = allprojs[n-1-k]（反转）

节点索引:
  用 (level, idx) 对标识节点，level 为层号，idx 为 levels[level] 内的下标
  不设虚根节点，levels[0] 直接为第 0 层节点
  cad_node 的 parent 字段为 idx（上一层内的下标），Level 0 节点无 parent
  cad_node 的 children 字段为 vector<size_t>（下一层内的下标列表）

单个 cell 的表示:
  vector<size_t> path，path[k] = levels[k] 中的下标
  从 path 可直接读取每层的半代数集描述和采样点
```

**功能规约**：

```
功能描述: 管理 CAD 分解的完整树形结构

前置条件:
  - _level_vars 非空
  - _level_polys[k] 中所有多项式含 _level_vars[k]（且不含 _level_vars[0..k-1] 的变量）

后置条件:
  - 提升完成后，每个叶节点路径对应 R^n 的一个胞腔
  - 每个叶节点可提取完整采样点 (α_1,...,α_n)

不变式:
  - _level_polys[k] 只增不减（可追加，已有索引不失效）
  - _levels[k] 的节点数 = 该层 Section 数 + Sector 数
  - _levels.size() == _level_vars.size() == _level_polys.size()
  - _level_vars[k]、_level_polys[k]、_levels[k] 对应同一个 CAD 层
```

### 2.2 `cad_node` — 树节点（内部类型）

```
结构名: cad_node

存储:
  type            : enum { SECTION, SECTOR }
  parent          : size_t              — 上一层内的下标（Level 0 节点无 parent）
  children        : vector<size_t>      — 下一层内的子节点下标列表

  // 根描述（Section 用 1 个，Sector 用 2 个）
  roots           : vector<cad_root>    — Section: size()=1, Sector: size()=2
  sample_point    : QQ                  — 有理采样点（Sector 必有；Section 近期无，远期代数数）

辅助类型:

  cad_root:
    is_inf        : int                 — -1 表示 -∞，+1 表示 +∞，0 表示有限根
    poly_idx      : size_t              — level_polys[level] 中的多项式索引（is_inf=0 时有效）
    root_idx      : size_t              — 该多项式的第几个实根（is_inf=0 时有效）
```

**功能规约**：

```
功能描述: 表示 CAD 树的单个节点（一维胞腔片段）

节点语义:
  - SECTION: roots 含 1 个 cad_root，标识一个代数根 {α}
             近期不存储采样点（Open CAD 不使用 Section）
             远期 sample_point 为三角列表示的代数数
  - SECTOR: roots 含 2 个 cad_root [left, right]，标识开区间 (left, right)
            left/right 可为 ±∞ 或有限根
            sample_point 为该开区间内的有理数

前置条件:
  - SECTION: roots.size() == 1, roots[0].is_inf == 0
  - SECTION: roots[0].poly_idx < level_polys[所在层].size()
  - SECTOR: roots.size() == 2
  - SECTOR: 有限根的 poly_idx/root_idx 有效
  - Level > 0 的节点: parent 为上一层的有效下标

构造保证（由提升过程维护，非静态可检查）:
  - SECTOR: 在父节点路径确定的采样点下，roots[0] 代表的实数 < roots[1] 代表的实数
  - SECTOR: sample_point 严格位于 roots[0] 和 roots[1] 对应的实数值之间

后置条件:
  - 叶节点（children 为空）代表一个完整胞腔
  - 从叶节点到根的路径上，完整描述该胞腔的半代数集表示
  - Sector 叶路径上的 sample_point 构成完整有理采样点向量
```

## 3. 接口规约

### 3.1 `cad_tree` 操作接口

#### 3.1.1 构造

```
cad_tree(vars, allprojs)
  输入: vars — 变量序 vector<variable>
        allprojs — __project_full 的输出 vector<vector<poly_type>>
  输出: 初始化的 cad_tree（仅含虚根节点，尚未提升）

  协议约定:
    - 调用方保证 vars 和 allprojs 一致（allprojs[k] 的多项式含 vars[k]）
    - cad_tree 接管 allprojs 的所有权（move 语义）
```

#### 3.1.2 节点操作（供算法函数调用）

```
add_sector(level, parent_idx, left_root, right_root, sample_point) → size_t
  功能: 在 levels[level] 添加一个 Sector 节点，parent_idx 为上一层内的下标
  返回: 新节点在 levels[level] 中的下标
  注: Level 0 节点无 parent，parent_idx 忽略

add_section(level, parent_idx, root) → size_t
  功能: 在 levels[level] 添加一个 Section 节点
  返回: 新节点在 levels[level] 中的下标

add_level_poly(level, poly) → size_t
  功能: 向 level_polys[level] 追加一个多项式
  返回: 新多项式在该层中的索引
```

#### 3.1.3 查询接口

```
get_sample_point(level, idx) → vector<QQ>
  从 (level, idx) 沿 parent 链回溯，收集各层采样点

node(level, idx) → const cad_node&
  访问 _levels[level][idx]

num_levels() → size_t
  返回层数 n

level_size(k) → size_t
  返回第 k 层的节点数

level_polys(k) → const vector<poly_type>&
  返回第 k 层的投影多项式集（提升序）

level_var(k) → const variable&
  返回第 k 层对应的变量（提升序）
```

注: `cad_root` 为公开类型，算法函数需要构造 `cad_root` 传给 `add_sector`/`add_section`。

### 3.2 算法函数接口

#### 3.2.1 __open_cad（字典序）

```
接口: __open_cad(polys, vars, method) → cad_tree

功能: 计算 Open CAD（仅 Sector 节点），要求输入多项式为 lex_<var_order> 字典序

输入:
  polys  — 输入多项式集 vector<polynomial_<ZZ, lex_<var_order>>>
  vars   — 变量序 vector<variable>，与 var_order 一致
  method — 投影算子选择（默认 LAZARD）

输出: 构建完成的 cad_tree，所有叶节点为 Sector

内部流程:
  1. 调用 __project_full(polys, vars, method) 得到 allprojs
  2. 构造 cad_tree(vars, allprojs)
  3. 从 Level 0 开始逐层提升:
     a. 取当前层投影多项式，代入父节点路径的采样点 → 一元多项式集
     b. 调用 realroot → 获取实根
     c. 在相邻根之间（含 ±∞ 边界）取有理采样点
     d. 调用 add_sector 添加子节点
  4. 返回 cad_tree
```

#### 3.2.2 open_cad（一般序，远期）

```
接口: open_cad(polys, vars, method) → cad_tree

功能: 一般序接口，先转换为字典序 lex_<custom_var_order>，再调用 __open_cad
状态: 远期实现，模式与现有 project() / __project() 一致
```

#### 3.2.3 代入求值（算法内部辅助）

```
接口: __eval_level_polys(tree, level, sample_path) → vector<upolynomial_ZZ>

功能: 将 tree.level_polys(level) 中的多项式，
      按 sample_path 代入前面各层变量的采样值，得到一元多项式集
输入: tree — cad_tree 引用
      level — 层号
      sample_path — 从根到当前节点的采样点序列
输出: 一元多项式向量（关于 vars[level] 的一元多项式）

协议约定:
  - sample_path.size() 等于需要代入的变量数
  - 返回的一元多项式可直接传给 realroot()
```

#### 3.2.4 __full_cad（远期）

```
接口: __full_cad(polys, vars, method) → cad_tree

功能: 计算完整 CAD（含 Section 和 Sector），字典序
状态: 远期实现，依赖代数数类型
```

## 4. 关键设计决策

### 4.1 分层存储

**决策**: 节点按层存储在 `vector<vector<cad_node>> levels` 中，`levels[k]` 为第 k 层所有节点。节点间通过层内下标引用。

**理由**:
- 提升算法逐层处理，分层存储下同层节点连续，缓存友好
- 单个 cell = `vector<size_t> path`（每层一个下标），表示轻量自然
- 部分 cell 集合（如加约束后的子集）= 每层一个有效下标集合
- `level_size(k)` 等按层查询 O(1)
- 不需要虚根节点，结构更简洁

### 4.2 采样点策略

**决策**: Sector 节点存有理采样点 `QQ`；Section 节点近期不存采样点，远期为三角列表示的代数数。

**理由**:
- Open CAD 只使用 Sector，有理采样点足够
- Section 的 `cad_root` 索引已完整标识代数根，远期实现代数数类型后按需计算
- 精确代数数表示（正规三角列 + 隔离区间列）留作远期扩展

### 4.3 level_polys 可追加

**决策**: `level_polys[k]` 用 `vector`，支持 `push_back`。

**理由**:
- 提升阶段处理 nullification 时需向当前层追加多项式
- 已有节点的 `poly_idx` 不会失效（vector 只追加）

### 4.4 数据结构与算法分离

**决策**: `cad_tree` 只做数据管理（存储节点、投影多项式、提供操作/查询接口），提升算法（`open_cad`、`full_cad`）为外部自由函数。

**理由**:
- CAD 有多种变体（Open/Full/Partial），都操作同一棵树，不应耦合进类
- 与 CLPoly 现有风格一致（`__project_full`、`realroot` 等均为自由函数）
- 可组合：投影和提升可独立调用，中途可追加多项式、切换策略

### 4.5 Open CAD 优先

**决策**: 首版只实现 Open CAD（仅 Sector 节点，不含 Section）。

**理由**:
- Open CAD 是最常用的应用场景（判断多项式系统在开胞腔上的符号）
- 实现简单，可快速验证整体流程正确性
- Full CAD（含 Section）作为后续扩展

## 5. 与现有代码的衔接

| 现有设施 | 在 cad_tree 中的使用 |
|---------|---------------------|
| `__project_full()` | 输出直接作为 `level_polys` 的初始值 |
| `realroot()` | 每层提升时调用，获取实根和隔离区间 |
| `sample_open_intervals()` | 实验版，未测试；实现 `__open_cad` 时根据需要修改或重写 |
| `uroot` | 提升过程中的中间结果，提取 `left()`/`right()` 构造采样点和 `cad_root` |
| `polynomial_<ZZ, lex_<var_order>>` | `level_polys` 的元素类型 |
| `assign()` (多项式求值) | `__eval_level_polys` 中代入采样值 |
