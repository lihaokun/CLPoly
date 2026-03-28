# CAD Cell 细化设计文档

> 基于 [architecture.md](architecture.md) 的架构设计，逐模块列出需要实现的函数。

## 模块 1: `cad_root` — 代数根引用

**文件**: `clpoly/cad_tree.hh`

```cpp
namespace clpoly {

class cad_root {
private:
    int _is_inf = 0;
    size_t _poly_idx = 0;
    size_t _root_idx = 0;    // 该多项式的第几个不同实根（不计重数，按值排序）
public:
    cad_root() {}
    cad_root(size_t poly_idx, size_t root_idx)
        : _is_inf(0), _poly_idx(poly_idx), _root_idx(root_idx) {}

    static cad_root inf()    { cad_root r; r._is_inf = 1; return r; }
    static cad_root neginf() { cad_root r; r._is_inf = -1; return r; }

    bool isinf() const    { return _is_inf == 1; }
    bool isneginf() const { return _is_inf == -1; }
    bool isfinite() const { return _is_inf == 0; }

    size_t poly_idx() const { return _poly_idx; }
    size_t root_idx() const { return _root_idx; }

    friend std::ostream& operator<<(std::ostream& stream, const cad_root& r);
};
```

---

## 模块 2: `cad_node` — 树节点

**文件**: `clpoly/cad_tree.hh`（与 cad_root 同文件）

```cpp
class cad_node {
private:
    enum cad_node_type { SECTION = 0, SECTOR = 1 };
    cad_node_type _type;
    size_t _parent = SIZE_MAX;          // 上一层内的下标，SIZE_MAX 表示无 parent（Level 0）
    size_t _child_begin = 0;            // 子节点在下一层中的起始下标
    size_t _child_count = 0;            // 子节点个数（子节点在下一层中连续存储）
    std::vector<cad_root> _roots;  // SECTION: size()=1, SECTOR: size()=2
    QQ _sample_point;
public:
    // SECTOR 构造
    cad_node(size_t parent, cad_root left, cad_root right, QQ sample)
        : _type(SECTOR), _parent(parent), _roots({left, right}), _sample_point(std::move(sample)) {}

    // SECTION 构造
    cad_node(size_t parent, cad_root root)
        : _type(SECTION), _parent(parent), _roots({root}) {}

    bool is_section() const { return _type == SECTION; }
    bool is_sector() const  { return _type == SECTOR; }

    size_t parent() const { return _parent; }
    size_t child_begin() const { return _child_begin; }
    size_t& child_begin() { return _child_begin; }
    size_t child_count() const { return _child_count; }
    size_t& child_count() { return _child_count; }
    const std::vector<cad_root>& roots() const { return _roots; }
    // todo: SECTION 的 sample_point 暂未定义，未来 Full CAD 需要处理代数数采样点
    const QQ& sample_point() const { assert(is_sector()); return _sample_point; }

    friend std::ostream& operator<<(std::ostream& stream, const cad_node& n);
};
```

---

## 模块 3: `cad_tree` — 管理类

**文件**: `clpoly/cad_tree.hh`

### 3.1 类定义

```cpp
template <class var_order>
class cad_tree {
public:
    using poly_type = polynomial_<ZZ, lex_<var_order>>;

private:
    std::vector<variable> _level_vars;    // 提升序变量
    std::vector<std::vector<poly_type>> _level_polys;     // 提升序投影多项式
    std::vector<std::vector<cad_node>> _levels;           // 分层节点

public:
    // 构造
    cad_tree() {}
    cad_tree(std::vector<variable> vars, std::vector<std::vector<poly_type>> allprojs);

    // 节点操作
    size_t add_sector(size_t level, size_t parent_idx, cad_root left, cad_root right, QQ sample);
    size_t add_section(size_t level, size_t parent_idx, cad_root root);
    size_t add_level_poly(size_t level, poly_type poly);  // 向指定层追加多项式

    // 查询
    const cad_node& node(size_t level, size_t idx) const { return _levels[level][idx]; }
    size_t num_levels() const { return _level_vars.size(); }       // 返回层数 n
    size_t level_size(size_t k) const { return _levels[k].size(); } // 第 k 层节点数
    const std::vector<poly_type>& level_polys(size_t k) const { return _level_polys[k]; } // 第 k 层投影多项式
    const variable& level_var(size_t k) const { return _level_vars[k]; } // 第 k 层对应的变量
    std::vector<QQ> get_sample_point(size_t level, size_t idx) const;    // 沿 parent 链回溯收集采样点

    // I/O
    template <class vo>
    friend std::ostream& operator<<(std::ostream& stream, const cad_tree<vo>& tree);
};
```

### 3.2 实现

```cpp
cad_tree(vars, allprojs):
    // 将投影阶段的字典序转为提升序（反转）
    n = vars.size()
    _level_vars.resize(n)
    _level_polys.resize(n)
    for k = 0 to n-1:
      _level_vars[k] = vars[n-1-k]
      _level_polys[k] = move(allprojs[n-1-k])
    _levels.resize(n)  // 预分配 n 层空 vector

add_sector(level, parent_idx, left, right, sample) → size_t:
    // 前置条件: level < num_levels()
    idx = _levels[level].size()
    _levels[level].emplace_back(parent_idx, left, right, std::move(sample))
    if (level > 0):  // Level 0 无 parent
      auto& p = _levels[level-1][parent_idx]
      if p.child_count() == 0:
        p.child_begin() = idx
      p.child_count()++
    return idx

add_section(level, parent_idx, root) → size_t:
    // 前置条件: level < num_levels()
    idx = _levels[level].size()
    _levels[level].emplace_back(parent_idx, root)
    if (level > 0):
      auto& p = _levels[level-1][parent_idx]
      if p.child_count() == 0:
        p.child_begin() = idx
      p.child_count()++
    return idx

add_level_poly(level, poly) → size_t:
    // 前置条件: level < num_levels()
    // 前置条件: poly 的主变量（first var）== _level_vars[level]
    // 追加多项式，已有索引不失效
    idx = _level_polys[level].size()
    _level_polys[level].push_back(move(poly))
    return idx

get_sample_point(level, idx) → vector<QQ>:
    // 从 (level, idx) 沿 parent 链回溯到 level 0
    path = vector<QQ>(level + 1)
    cur_idx = idx
    for k = level downto 0:
      path[k] = _levels[k][cur_idx].sample_point()
      if k > 0:
        cur_idx = _levels[k][cur_idx].parent()
    return path
```

---

## 模块 4: `__open_cad` — Open CAD 算法

**文件**: `clpoly/cad.hh`（添加到现有文件）

### 4.1 主函数

```cpp
template <class var_order>
cad_tree<var_order> __open_cad(
    const vector<polynomial_<ZZ, lex_<var_order>>>& polys,
    const vector<variable>& vars,
    projection_method method = projection_method::LAZARD
)
  功能: 计算 Open CAD（仅 Sector）

  实现:
    // 1. 投影
    allprojs = __project_full(polys, vars, method)

    // 2. 构造树（构造函数内部反转为提升序）
    tree = cad_tree<var_order>(vars, move(allprojs))
    n = vars.size()

    // 3. Level 0: _level_polys[0] 为一元多项式，无需代入
    auto [roots_0, info_0] = realroot(tree.level_polys(0))
    __lift_open_level(tree, 0, SIZE_MAX, roots_0, info_0)

    // 4. Level 1 到 Level n-1: 逐层提升
    for level = 1 to n-1:
      for parent_idx = 0 to tree.level_size(level-1)-1:
        eval_polys = __eval_level_polys(tree, level, tree.get_sample_point(level-1, parent_idx))
        auto [roots_k, info_k] = realroot(eval_polys)
        __lift_open_level(tree, level, parent_idx, roots_k, info_k)

    return tree
```

### 4.2 辅助函数: 单层提升

```cpp
template <class var_order>
void __lift_open_level(
    cad_tree<var_order>& tree,
    size_t level,
    size_t parent_idx,
    const vector<uroot>& roots,
    const vector<vector<pair<uint64_t,uint64_t>>>& poly_mult_info
)
  功能: 在指定 level 的 parent 下根据实根创建所有 Sector 节点
        Level 0 时 parent_idx = SIZE_MAX（无 parent）

  实现:
    // 根数为 0: 整条实数线是一个 Sector (-∞, +∞)
    if roots.empty():
      tree.add_sector(level, parent_idx, cad_root::neginf(), cad_root::inf(), QQ(0))
      return

    // 预计算所有根对应的 cad_root
    cad_roots = __make_cad_roots(poly_mult_info, tree.level_polys(level).size())

    // (-∞, roots[0])
    tree.add_sector(level, parent_idx, cad_root::neginf(), cad_roots[0],
                    roots[0].left() - 1)

    // (roots[i], roots[i+1]) for i = 0..m-2
    for i = 0 to roots.size()-2:
      assert(roots[i].right() < roots[i+1].left())  // 隔离区间不重叠
      tree.add_sector(level, parent_idx, cad_roots[i], cad_roots[i+1],
                      (roots[i].right() + roots[i+1].left()) / 2)

    // (roots[m-1], +∞)
    tree.add_sector(level, parent_idx, cad_roots.back(), cad_root::inf(),
                    roots.back().right() + 1)
```

### 4.3 辅助函数: 代入求值

```cpp
template <class var_order>
vector<polynomial_<ZZ, lex_<var_order>>> __eval_level_polys(
    const cad_tree<var_order>& tree,
    size_t level,
    const vector<QQ>& sample_path
)
  功能: 将 _level_polys[level] 代入采样点，得到低维 ZZ 多项式
  实现:
    // 构造代入映射: sample_path[k] 对应 level k 的变量 _level_vars[k]
    // sample_path 包含 level 0 到 level-1 的采样点（之前层已确定的值）
    subst_map = map<variable, QQ>{}
    for k = 0 to sample_path.size()-1:
      subst_map[tree.level_var(k)] = sample_path[k]

    result = {}
    for each poly in tree.level_polys(level):
      // 一次性代入所有已确定变量: ZZ 多项式 → QQ 多项式
      p_qq = assign<QQ, ZZ, QQ>(poly, subst_map)

      // QQ 多项式转 ZZ 多项式（乘公分母，不影响实根）
      poly_convert(p_qq, p_zz)
      result.push_back(p_zz)

    return result

  复用: polynomial.hh 的 assign<Tc1,Tc2,Tc3,To>(poly, map)
        polynomial_convert.hh 的 poly_convert(polynomial_<QQ> → polynomial_<ZZ>)
```

### 4.4 辅助函数: realroot 结果到 cad_root 映射

```cpp
vector<cad_root> __make_cad_roots(
    const vector<vector<pair<uint64_t,uint64_t>>>& poly_mult_info,
    size_t num_polys
)
  功能: 一次性为所有根构造 cad_root，O(n)
  实现:
    // poly_root_count[p] = 多项式 p 已出现的不同实根数（不计重数）
    poly_root_count = vector<size_t>(num_polys, 0)
    result = vector<cad_root>(poly_mult_info.size())

    for i = 0 to poly_mult_info.size()-1:
      // 选第一个多项式作为 cad_root 的标识
      poly_idx = poly_mult_info[i][0].first
      local_idx = poly_root_count[poly_idx]
      result[i] = cad_root(poly_idx, local_idx)

      // 对该根所属的所有多项式递增计数，保证后续根的局部索引正确
      for each (pidx, _) in poly_mult_info[i]:
        poly_root_count[pidx]++

    return result
```

---

## 模块 5: 测试

**文件**: `test/test_cad_tree.cc`, `test/test_open_cad.cc`

### 5.1 基本功能测试

```
test_cad_root_creation
  测试 cad_root::inf(), cad_root::neginf(), cad_root(poly_idx, root_idx) 的构造

test_cad_tree_construction
  构造空 cad_tree，检查 num_levels, level_size

test_cad_tree_add_sector
  手动添加节点，检查 parent/children 关系、node 查询

test_cad_tree_get_sample_point
  多层手动构建树，验证 get_sample_point 的回溯正确性
```

### 5.2 集成测试: __open_cad

```
test_open_cad_univariate
  单变量多项式（n=1），验证 Sector 节点数 = 实根数 + 1

test_open_cad_bivariate
  二元多项式（如 x^2 + y^2 - 1），验证树结构和采样点

test_open_cad_trivariate
  三元多项式，验证完整 3 层树

test_open_cad_no_real_roots
  无实根的多项式，验证只有一个 Sector (-∞, +∞)
```

---

## 与已有代码的复用点

| 函数 | 来源 | 用途 |
|------|------|------|
| `__project_full` | cad.hh | __open_cad 第一步 |
| `realroot` | realroot.hh | 每层提升获取实根 |
| `assign<QQ,ZZ,QQ>` | polynomial.hh | __eval_level_polys 代入 QQ 采样值 |
| `poly_convert(QQ→ZZ)` | polynomial_convert.hh | QQ 多项式乘公分母转回 ZZ |
| `uroot` 比较/区间 | realroot.hh | 采样点计算 |

## 待确认问题

（已全部解决）
