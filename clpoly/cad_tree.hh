/**
 * @file cad_tree.hh
 * @author ntimesp(nxp@mail.ustc.edu.cn)
 * @brief CAD 树形数据结构：cad_root, cad_node, cad_tree
 *
 * cad_root — 代数根的索引引用，标识某层某多项式的第几个不同实根（或 ±∞）。
 *            纯值类型，不持有实际代数数值，仅用于标注 cad_node 的边界。
 * cad_node — CAD 树的节点，表示一维胞腔片段（Sector 或 Section）。
 *            记录 parent 索引、children 范围、边界根引用和采样点。
 * cad_tree — 整棵 CAD 分解树的容器，按层（level）存储节点。
 *            持有提升序的变量列表、各层投影多项式和节点数组。
 *            目前构建完成后为只读结构（children 用 child_begin + child_count
 *            连续存储，要求同一 parent 的子节点按顺序插入）。
 */
#ifndef CLPOLY_CAD_TREE_HH
#define CLPOLY_CAD_TREE_HH

#include <clpoly/polynomial.hh>
#include <clpoly/number.hh>
#include <vector>
#include <iostream>
#include <climits>
#include <cassert>

namespace clpoly{

    // ========== cad_root ==========
    // 代数根的索引引用：标识 level_polys[level][poly_idx] 的第 root_idx 个不同实根，
    // 或表示 ±∞。

    class cad_root {
    private:
        int _is_inf = 0;       // -1: -∞, 0: 有限根, +1: +∞
        size_t _poly_idx = 0;  // level_polys 中的多项式下标（_is_inf==0 时有效）
        size_t _root_idx = 0;  // 该多项式的第几个不同实根（不计重数，按值排序）
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

        friend bool operator==(const cad_root& a, const cad_root& b)
        {
            if (a._is_inf != b._is_inf) return false;
            if (a._is_inf != 0) return true;  // 都是同符号的 ∞
            return a._poly_idx == b._poly_idx && a._root_idx == b._root_idx;
        }
        friend bool operator!=(const cad_root& a, const cad_root& b) { return !(a == b); }

        friend std::ostream& operator<<(std::ostream& stream, const cad_root& r)
        {
            if (r.isinf())
            {
                stream << "+inf";
                return stream;
            }
            if (r.isneginf())
            {
                stream << "-inf";
                return stream;
            }
            stream << "(poly=" << r._poly_idx << ",root=" << r._root_idx << ")";
            return stream;
        }
    };

    // ========== cad_node ==========
    // CAD 树的节点，表示一维胞腔片段。
    // SECTION: 位于一个代数根上（roots 含 1 个 cad_root）
    // SECTOR:  位于两个相邻根之间的开区间（roots 含 2 个 cad_root，可为 ±∞）

    class cad_node {
    private:
        enum cad_node_type { SECTION = 0, SECTOR = 1 };
        cad_node_type _type;
        size_t _parent = SIZE_MAX;          // 上一层内的下标，SIZE_MAX 表示无 parent（Level 0）
        size_t _child_begin = 0;             // 子节点在下一层中的起始下标
        size_t _child_count = 0;             // 子节点个数（子节点在下一层中连续存储）
        std::vector<cad_root> _roots;       // SECTION: size()=1, SECTOR: size()=2
        QQ _sample_point;                   // 仅 SECTOR 使用
    public:
        // SECTOR 构造: 开区间 (left, right) 内取 sample 为采样点
        cad_node(size_t parent, cad_root left, cad_root right, QQ sample)
            : _type(SECTOR), _parent(parent), _roots({left, right}), _sample_point(std::move(sample)) {}

        // SECTION 构造: 位于 root 标识的代数根上
        cad_node(size_t parent, cad_root root)
            : _type(SECTION), _parent(parent), _roots({root})
        {
            assert(root.isfinite());  // Section 的根不能是 ±∞
        }

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

        friend bool operator==(const cad_node& a, const cad_node& b)
        {
            if (a._type != b._type) return false;
            if (a._parent != b._parent) return false;
            if (a._child_begin != b._child_begin) return false;
            if (a._child_count != b._child_count) return false;
            if (a._roots.size() != b._roots.size()) return false;
            for (size_t i = 0; i < a._roots.size(); ++i)
                if (a._roots[i] != b._roots[i]) return false;
            if (a._type == SECTOR && a._sample_point != b._sample_point) return false;
            return true;
        }
        friend bool operator!=(const cad_node& a, const cad_node& b) { return !(a == b); }

        friend std::ostream& operator<<(std::ostream& stream, const cad_node& n)
        {
            if (n.is_section())
            {
                stream << "Section{" << n._roots[0] << "}";
            }
            else
            {
                stream << "Sector{" << n._roots[0] << ", " << n._roots[1]
                       << ", sample=" << n._sample_point << "}";
            }
            return stream;
        }
    };

    // ========== cad_tree ==========
    // CAD 分解的树形管理结构。
    //
    // 内部使用提升序索引：
    //   level 0 = 一元多项式层（最先提升）
    //   level n-1 = 最后提升（叶层）
    //   _level_vars[k]、_level_polys[k]、_levels[k] 三者对应同一层
    //
    // 外部输入（vars, allprojs）按字典序，构造函数内部反转为提升序。

    template <class var_order>
    class cad_tree {
    public:
        using poly_type = polynomial_<ZZ, lex_<var_order>>;

    private:
        std::vector<variable> _level_vars;                // 提升序变量
        std::vector<std::vector<poly_type>> _level_polys; // 提升序投影多项式（可追加）
        std::vector<std::vector<cad_node>> _levels;       // 分层节点

    public:
        // 默认构造，空树
        cad_tree() {}

        // 从投影结果构造，内部反转为提升序
        // vars: 字典序变量（vars[0] = lex 最小）
        // allprojs: __project_full 的输出（allprojs[k] 含 vars[k]）
        cad_tree(std::vector<variable> vars, std::vector<std::vector<poly_type>> allprojs)
        {
            assert(vars.size() == allprojs.size());
            size_t n = vars.size();
            _level_vars.resize(n);
            _level_polys.resize(n);
            for (size_t k = 0; k < n; ++k)
            {
                _level_vars[k] = vars[n-1-k];
                _level_polys[k] = std::move(allprojs[n-1-k]);
            }
            _levels.resize(n);
        }

        // ---- 节点操作 ----

        // 在 level 层添加一个 SECTOR 节点
        // Level 0 时 parent_idx 应为 SIZE_MAX（无 parent）
        size_t add_sector(size_t level, size_t parent_idx, cad_root left, cad_root right, QQ sample)
        {
            assert(level < num_levels());
            assert(level == 0 || parent_idx < _levels[level-1].size());
            size_t idx = _levels[level].size();
            _levels[level].emplace_back(parent_idx, left, right, std::move(sample));
            if (level > 0)
            {
                auto& p = _levels[level-1][parent_idx];
                if (p.child_count() == 0)
                    p.child_begin() = idx;
                p.child_count()++;
            }
            return idx;
        }

        // 在 level 层添加一个 SECTION 节点
        size_t add_section(size_t level, size_t parent_idx, cad_root root)
        {
            assert(level < num_levels());
            assert(level == 0 || parent_idx < _levels[level-1].size());
            size_t idx = _levels[level].size();
            _levels[level].emplace_back(parent_idx, root);
            if (level > 0)
            {
                auto& p = _levels[level-1][parent_idx];
                if (p.child_count() == 0)
                    p.child_begin() = idx;
                p.child_count()++;
            }
            return idx;
        }

        // 向指定层追加一个投影多项式
        size_t add_level_poly(size_t level, poly_type poly)
        {
            assert(level < num_levels());
            size_t idx = _level_polys[level].size();
            _level_polys[level].push_back(std::move(poly));
            return idx;
        }

        // ---- 查询 ----

        const cad_node& node(size_t level, size_t idx) const
        {
            assert(level < num_levels() && idx < _levels[level].size());
            return _levels[level][idx];
        }
        size_t num_levels() const { return _level_vars.size(); }
        size_t level_size(size_t k) const
        {
            assert(k < num_levels());
            return _levels[k].size();
        }
        const std::vector<poly_type>& level_polys(size_t k) const
        {
            assert(k < num_levels());
            return _level_polys[k];
        }
        const variable& level_var(size_t k) const
        {
            assert(k < num_levels());
            return _level_vars[k];
        }

        // 沿 parent 链回溯，收集从 level 0 到 level 的各层采样点
        std::vector<QQ> get_sample_point(size_t level, size_t idx) const
        {
            assert(level < num_levels() && idx < _levels[level].size());
            std::vector<QQ> path(level + 1);
            size_t cur_idx = idx;
            for (size_t k = level + 1; k > 0; --k)
            {
                path[k-1] = _levels[k-1][cur_idx].sample_point();
                if (k-1 > 0)
                    cur_idx = _levels[k-1][cur_idx].parent();
            }
            return path;
        }

        // I/O
        template <class vo>
        friend std::ostream& operator<<(std::ostream& stream, const cad_tree<vo>& tree);

        // 比较
        friend bool operator==(const cad_tree& a, const cad_tree& b)
        {
            if (a._level_vars != b._level_vars) return false;
            if (a._level_polys != b._level_polys) return false;
            if (a._levels.size() != b._levels.size()) return false;
            for (size_t k = 0; k < a._levels.size(); ++k)
            {
                if (a._levels[k].size() != b._levels[k].size()) return false;
                for (size_t i = 0; i < a._levels[k].size(); ++i)
                    if (a._levels[k][i] != b._levels[k][i]) return false;
            }
            return true;
        }
        friend bool operator!=(const cad_tree& a, const cad_tree& b) { return !(a == b); }
    };

    template <class var_order>
    std::ostream& operator<<(std::ostream& stream, const cad_tree<var_order>& tree)
    {
        stream << "cad_tree: " << tree.num_levels() << " levels\n";
        for (size_t k = 0; k < tree.num_levels(); ++k)
        {
            stream << "  Level " << k << " (var=" << tree._level_vars[k]
                   << ", polys=" << tree._level_polys[k].size()
                   << ", nodes=" << tree._levels[k].size() << "):\n";
            for (size_t i = 0; i < tree._levels[k].size(); ++i)
            {
                const auto& nd = tree._levels[k][i];
                stream << "    [" << i << "] " << nd;
                // parent 信息
                if (nd.parent() != SIZE_MAX)
                    stream << " parent=" << nd.parent();
                // children 信息
                if (nd.child_count() > 0)
                {
                    stream << " children={";
                    for (size_t j = 0; j < nd.child_count(); ++j)
                    {
                        if (j > 0) stream << ",";
                        stream << nd.child_begin() + j;
                    }
                    stream << "}";
                }
                stream << "\n";
            }
        }
        return stream;
    }

} // namespace clpoly
#endif
