/**
 * @file polynomial_factorize_univar.hh
 * @brief ZZ[x] 单变量因式分解：ZZ 辅助函数、Hensel 提升 (M2)、
 *        Zassenhaus 重组 (M3)、素数选择、factorize(upolynomial_<ZZ>)。
 *
 * 本文件可独立包含。供 polynomial_factorize_wang.hh 和
 * polynomial_factorize.hh 使用。
 *
 * 注意：factorization<> 结构体定义在本文件中。
 */
#ifndef CLPOLY_POLYNOMIAL_FACTORIZE_UNIVAR_HH
#define CLPOLY_POLYNOMIAL_FACTORIZE_UNIVAR_HH

#include <clpoly/polynomial_factorize_zp.hh>
#include <clpoly/polynomial.hh>
#include <clpoly/polynomial_gcd.hh>
#include <boost/math/special_functions/prime.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>

namespace clpoly{

    // ================================================================
    // §8.1 factorization 返回类型
    // ================================================================
    template<class Poly>
    struct factorization {
        typename Poly::coeff_type content;
        std::vector<std::pair<Poly, uint64_t>> factors;
    };

    // ================================================================
    // §4.2 对称模约化
    // ================================================================

    // §4.2.1 标量对称模: 将 a 约化到 (-m/2, m/2]
    inline ZZ __symmetric_mod(const ZZ& a, const ZZ& m)
    {
        ZZ r;
        ZZ::fdiv_r(r, a, m);   // r ∈ [0, m)
        ZZ half;
        ZZ::fdiv_q(half, m, ZZ(2));
        if (r > half)
            r -= m;
        return r;
    }

    // §4.2.2 多项式对称模: 对每个系数做对称模约化
    inline upolynomial_<ZZ> __upoly_symmetric_mod(
        const upolynomial_<ZZ>& f,
        const ZZ& m)
    {
        upolynomial_<ZZ> result;
        result.reserve(f.size());
        for (auto& term : f)
        {
            ZZ c = __symmetric_mod(term.second, m);
            if (c)
                result.push_back(std::make_pair(term.first, std::move(c)));
        }
        return result;
    }

    // ================================================================
    // §4.3 范数计算
    // ================================================================

    // L2 范数的平方 (系数平方和)
    inline ZZ __upoly_norm_l2_sq(const upolynomial_<ZZ>& f)
    {
        ZZ s(0);
        for (auto& term : f)
            s += term.second * term.second;
        return s;
    }

    // ================================================================
    // §4.4 Mignotte 界
    // ================================================================

    // 二项式系数 C(n, k) (ZZ 精确计算)
    inline ZZ __binomial(int64_t n, int64_t k)
    {
        if (k < 0 || k > n) return ZZ(0);
        if (k == 0 || k == n) return ZZ(1);
        if (k > n - k) k = n - k;
        ZZ result(1);
        for (int64_t i = 0; i < k; ++i)
        {
            result *= ZZ(n - i);
            result /= ZZ(i + 1);
        }
        return result;
    }

    // 整数平方根的上界 (ceiling of sqrt)
    inline ZZ __isqrt_ceil(const ZZ& n)
    {
        if (n <= ZZ(0)) return ZZ(0);
        // use sizeinbase to get initial estimate
        size_t bits = n.sizeinbase(2);
        ZZ x(1);
        x <<= (bits + 1) / 2;
        // Newton's method
        while (true)
        {
            ZZ x1 = (x + n / x) / ZZ(2);
            if (x1 >= x)
                break;
            x = x1;
        }
        // ensure x*x >= n
        if (x * x < n)
            x += 1;
        return x;
    }

    // Mignotte bound: B = C(n, ⌊n/2⌋) · ‖f‖₂
    inline ZZ __mignotte_bound(const upolynomial_<ZZ>& f)
    {
        assert(!f.empty());
        int64_t n = get_deg(f);
        ZZ binom = __binomial(n, n / 2);
        ZZ norm_sq = __upoly_norm_l2_sq(f);
        ZZ norm = __isqrt_ceil(norm_sq);
        return binom * norm;
    }

    // ================================================================
    // §4.1.4 ZZ 系数的模除法
    // ================================================================

    // 将 ZZ 多项式各系数 mod m (非负余数 [0, m))
    inline void __upoly_mod_coeff(upolynomial_<ZZ>& f, const ZZ& m)
    {
        auto it = f.data().begin();
        auto out = it;
        for (; it != f.data().end(); ++it)
        {
            ZZ::fdiv_r(it->second, it->second, m);
            if (it->second)
            {
                if (out != it)
                    *out = std::move(*it);
                ++out;
            }
        }
        f.data().erase(out, f.data().end());
    }

    // 计算 f = q·g + r in Z_m[x] (长除法, 每步 mod m)
    inline void __upoly_divmod_mod(
        upolynomial_<ZZ>& q,
        upolynomial_<ZZ>& r,
        const upolynomial_<ZZ>& f,
        const upolynomial_<ZZ>& g,
        const ZZ& m)
    {
        assert(!g.empty());
        q.clear();
        r = f;
        __upoly_mod_coeff(r, m);

        ZZ lc_inv;
        bool ok = ZZ::invert(lc_inv, g.front().second, m);
        assert(ok); (void)ok;

        int64_t deg_g = get_deg(g);

        while (!r.empty() && get_deg(r) >= deg_g)
        {
            int64_t d = get_deg(r) - deg_g;
            ZZ coeff = r.front().second * lc_inv;
            ZZ::fdiv_r(coeff, coeff, m);

            if (!coeff)
            {
                // leading term vanishes mod m, remove it
                r.data().erase(r.data().begin());
                continue;
            }

            q.push_back(std::make_pair(umonomial(d), coeff));

            // r -= coeff * x^d * g
            upolynomial_<ZZ> new_r;
            new_r.reserve(r.size());
            auto r_it = r.data().begin();
            auto g_it = g.begin();

            // skip r's leading term (it gets cancelled)
            ++r_it;
            ++g_it;

            while (r_it != r.data().end() || g_it != g.end())
            {
                if (g_it == g.end())
                {
                    ZZ c = r_it->second;
                    ZZ::fdiv_r(c, c, m);
                    if (c)
                        new_r.push_back(std::make_pair(r_it->first, std::move(c)));
                    ++r_it;
                }
                else if (r_it == r.data().end())
                {
                    int64_t deg_term = g_it->first.deg() + d;
                    ZZ c = m - (coeff * g_it->second % m) % m;
                    ZZ::fdiv_r(c, c, m);
                    if (c)
                        new_r.push_back(std::make_pair(umonomial(deg_term), std::move(c)));
                    ++g_it;
                }
                else
                {
                    int64_t deg_r = r_it->first.deg();
                    int64_t deg_term = g_it->first.deg() + d;
                    if (deg_r > deg_term)
                    {
                        ZZ c = r_it->second;
                        ZZ::fdiv_r(c, c, m);
                        if (c)
                            new_r.push_back(std::make_pair(r_it->first, std::move(c)));
                        ++r_it;
                    }
                    else if (deg_r < deg_term)
                    {
                        ZZ c = m - (coeff * g_it->second % m) % m;
                        ZZ::fdiv_r(c, c, m);
                        if (c)
                            new_r.push_back(std::make_pair(umonomial(deg_term), std::move(c)));
                        ++g_it;
                    }
                    else
                    {
                        ZZ c = r_it->second - coeff * g_it->second;
                        ZZ::fdiv_r(c, c, m);
                        if (c)
                            new_r.push_back(std::make_pair(r_it->first, std::move(c)));
                        ++r_it;
                        ++g_it;
                    }
                }
            }
            r.data() = std::move(new_r.data());
        }
    }

    // ================================================================
    // §4.5 Zp→ZZ 多项式转换
    // ================================================================

    inline upolynomial_<ZZ> __upoly_Zp_to_ZZ(const upolynomial_<Zp>& f)
    {
        upolynomial_<ZZ> result;
        result.reserve(f.size());
        for (auto& term : f)
            result.push_back(std::make_pair(term.first, ZZ(static_cast<int64_t>(term.second.number()))));
        return result;
    }

    // ================================================================
    // §6. M2: Hensel 提升
    // ================================================================

    // §6.3 Hensel 树节点
    struct __hensel_node {
        upolynomial_<ZZ> g;     // 左子树因子之积
        upolynomial_<ZZ> h;     // 右子树因子之积
        upolynomial_<ZZ> s;     // Bézout 系数: s·g + t·h ≡ 1 (mod m)
        upolynomial_<ZZ> t;     // Bézout 系数
        int left;               // 左子节点索引 (-1 = 叶子)
        int right;              // 右子节点索引 (-1 = 叶子)
        int leaf_start;         // 叶子范围起始
        int leaf_end;           // 叶子范围结束
    };

    // ZZ 多项式乘法并 mod m
    inline upolynomial_<ZZ> __upoly_mul_mod(
        const upolynomial_<ZZ>& a,
        const upolynomial_<ZZ>& b,
        const ZZ& m)
    {
        upolynomial_<ZZ> result = a * b;
        __upoly_mod_coeff(result, m);
        return result;
    }

    // §6.4 构建初始提升树
    inline void __hensel_tree_build_recursive(
        std::vector<__hensel_node>& nodes,
        const std::vector<upolynomial_<Zp>>& factors,
        uint32_t p,
        int start, int end,   // factors[start..end)
        int parent_idx)       // 在 nodes 中的索引
    {
        int mid = (start + end) / 2;

        // 计算 g = product of factors[start..mid)
        upolynomial_<Zp> g_zp = factors[start];
        for (int i = start + 1; i < mid; ++i)
            g_zp = g_zp * factors[i];

        // 计算 h = product of factors[mid..end)
        upolynomial_<Zp> h_zp = factors[mid];
        for (int i = mid + 1; i < end; ++i)
            h_zp = h_zp * factors[i];

        // 扩展 GCD: s·g + t·h = gcd mod p
        upolynomial_<Zp> s_zp, t_zp;
        polynomial_GCD(g_zp, h_zp, s_zp, t_zp);

        // Zp → ZZ
        nodes[parent_idx].g = __upoly_Zp_to_ZZ(g_zp);
        nodes[parent_idx].h = __upoly_Zp_to_ZZ(h_zp);
        nodes[parent_idx].s = __upoly_Zp_to_ZZ(s_zp);
        nodes[parent_idx].t = __upoly_Zp_to_ZZ(t_zp);
        nodes[parent_idx].leaf_start = start;
        nodes[parent_idx].leaf_end = end;

        // 递归左子树
        if (mid - start >= 2)
        {
            int left_idx = nodes.size();
            nodes.push_back({});
            nodes[parent_idx].left = left_idx;
            __hensel_tree_build_recursive(nodes, factors, p, start, mid, left_idx);
        }
        else
        {
            nodes[parent_idx].left = -1; // 叶子
        }

        // 递归右子树
        if (end - mid >= 2)
        {
            int right_idx = nodes.size();
            nodes.push_back({});
            nodes[parent_idx].right = right_idx;
            __hensel_tree_build_recursive(nodes, factors, p, mid, end, right_idx);
        }
        else
        {
            nodes[parent_idx].right = -1; // 叶子
        }
    }

    inline std::vector<__hensel_node> __hensel_tree_build(
        const std::vector<upolynomial_<Zp>>& factors,
        uint32_t p)
    {
        assert(factors.size() >= 2);
        std::vector<__hensel_node> nodes;
        nodes.push_back({}); // root at index 0
        __hensel_tree_build_recursive(nodes, factors, p, 0, (int)factors.size(), 0);
        return nodes;
    }

    // §6.5 单步二次 Hensel 提升
    inline void __hensel_step(
        __hensel_node& node,
        const upolynomial_<ZZ>& f,
        const ZZ& m)
    {
        ZZ m2 = m * m;

        // --- 第一部分: 提升因子 ---

        // 1. e = (f - g*h) / m, then e mod m
        upolynomial_<ZZ> gh = node.g * node.h;
        upolynomial_<ZZ> e = f - gh;
        // e 的每个系数精确整除 m
        for (auto& term : e.data())
        {
            ZZ::fdiv_q(term.second, term.second, m);
            ZZ::fdiv_r(term.second, term.second, m);
        }
        // 去掉零项
        auto it = e.data().begin();
        auto out = it;
        for (; it != e.data().end(); ++it)
        {
            if (it->second) { if (out != it) *out = std::move(*it); ++out; }
        }
        e.data().erase(out, e.data().end());

        // 2. se = s * e, divmod by h (mod m)
        upolynomial_<ZZ> se = node.s * e;
        upolynomial_<ZZ> q_se, r_se;
        __upoly_divmod_mod(q_se, r_se, se, node.h, m);

        // 3. 更新 g, h
        // tau = t*e + q_se*g, all mod m
        upolynomial_<ZZ> te = node.t * e;
        upolynomial_<ZZ> qg = q_se * node.g;
        upolynomial_<ZZ> tau = te + qg;
        __upoly_mod_coeff(tau, m);

        // g_new = g + m * tau (mod m²)
        for (auto& term : tau.data())
            term.second *= m;
        node.g = node.g + tau;
        __upoly_mod_coeff(node.g, m2);

        // h_new = h + m * r_se (mod m²)
        for (auto& term : r_se.data())
            term.second *= m;
        node.h = node.h + r_se;
        __upoly_mod_coeff(node.h, m2);

        // --- 第二部分: 提升 Bézout 系数 ---

        // 4. e' = (1 - s*g - t*h) / m, then e' mod m
        upolynomial_<ZZ> sg = node.s * node.g;
        upolynomial_<ZZ> th = node.t * node.h;
        // construct 1
        upolynomial_<ZZ> one;
        one.push_back(std::make_pair(umonomial(0), ZZ(1)));
        upolynomial_<ZZ> ep = one - sg - th;
        for (auto& term : ep.data())
        {
            ZZ::fdiv_q(term.second, term.second, m);
            ZZ::fdiv_r(term.second, term.second, m);
        }
        it = ep.data().begin();
        out = it;
        for (; it != ep.data().end(); ++it)
        {
            if (it->second) { if (out != it) *out = std::move(*it); ++out; }
        }
        ep.data().erase(out, ep.data().end());

        // 5. divmod(s*e', h, m)
        upolynomial_<ZZ> sep = node.s * ep;
        upolynomial_<ZZ> q_sep, r_sep;
        __upoly_divmod_mod(q_sep, r_sep, sep, node.h, m);

        // s_new = s + m * r_sep (mod m²)
        for (auto& term : r_sep.data())
            term.second *= m;
        node.s = node.s + r_sep;
        __upoly_mod_coeff(node.s, m2);

        // t_new = t + m * (t*e' + q_sep*g) (mod m²)
        upolynomial_<ZZ> tep = node.t * ep;
        upolynomial_<ZZ> qpg = q_sep * node.g;
        upolynomial_<ZZ> tau2 = tep + qpg;
        __upoly_mod_coeff(tau2, m);
        for (auto& term : tau2.data())
            term.second *= m;
        node.t = node.t + tau2;
        __upoly_mod_coeff(node.t, m2);
    }

    // 从 Hensel 树中提取叶子因子
    inline void __hensel_extract_factors(
        const std::vector<__hensel_node>& nodes,
        int idx,
        std::vector<upolynomial_<ZZ>>& factors)
    {
        const auto& node = nodes[idx];
        if (node.left == -1)
            factors.push_back(node.g);
        else
            __hensel_extract_factors(nodes, node.left, factors);

        if (node.right == -1)
            factors.push_back(node.h);
        else
            __hensel_extract_factors(nodes, node.right, factors);
    }

    // 递归自顶向下 Hensel 提升:
    // 先提升当前节点 (target ≡ g*h mod m → mod m²)
    // 再用更新后的 g, h 作为子节点的 target 递归提升
    inline void __hensel_lift_recursive(
        std::vector<__hensel_node>& nodes,
        int idx,
        const upolynomial_<ZZ>& target,
        const ZZ& m)
    {
        __hensel_step(nodes[idx], target, m);

        if (nodes[idx].left != -1)
            __hensel_lift_recursive(nodes, nodes[idx].left, nodes[idx].g, m);
        if (nodes[idx].right != -1)
            __hensel_lift_recursive(nodes, nodes[idx].right, nodes[idx].h, m);
    }

    // §6.6 多因子 Hensel 提升 (M2 入口)
    inline std::pair<std::vector<upolynomial_<ZZ>>, ZZ>
    __hensel_lift(
        const upolynomial_<ZZ>& f,
        const std::vector<upolynomial_<Zp>>& factors,
        uint32_t p)
    {
        assert(factors.size() >= 2);
        assert(!f.empty());

        // 1. 确定提升精度
        ZZ B = __mignotte_bound(f);
        ZZ lc_f = f.front().second;
        if (lc_f < ZZ(0)) lc_f = -lc_f;
        ZZ target = ZZ(2) * lc_f * B;  // 需要 p^k > target

        // 2. 处理首项系数: 将 lc(f) 分配到 factor[0]
        //    修改 factor[0]: 乘以 lc(f) mod p
        std::vector<upolynomial_<Zp>> factors_adj = factors;
        Zp lc_mod_p(f.front().second, p);
        for (auto& term : factors_adj[0])
            term.second *= lc_mod_p;
        factors_adj[0].normalization();

        // 3. 构建初始树 (使用调整后的因子)
        auto nodes = __hensel_tree_build(factors_adj, p);

        // 4. 二次提升 (自顶向下递归)
        ZZ m(p);
        while (m <= target)
        {
            __hensel_lift_recursive(nodes, 0, f, m);
            m = m * m;
        }

        // 5. 提取最终因子
        std::vector<upolynomial_<ZZ>> result;
        __hensel_extract_factors(nodes, 0, result);

        return {std::move(result), std::move(m)};
    }

    // ================================================================
    // §7. M3: 因子重组 (Zassenhaus)
    // ================================================================

    // §4.3 L1 范数 (系数绝对值之和)
    inline ZZ __upoly_norm_l1(const upolynomial_<ZZ>& f)
    {
        ZZ s(0);
        for (auto& term : f)
        {
            ZZ a = term.second;
            if (a < ZZ(0)) a = -a;
            s += a;
        }
        return s;
    }

    // §4.5 本原化: 提取内容并本原化
    inline std::pair<ZZ, upolynomial_<ZZ>> __upoly_primitive(upolynomial_<ZZ> f)
    {
        if (f.empty()) return {ZZ(1), std::move(f)};
        ZZ c = cont(f);
        if (f.front().second < ZZ(0)) c = -c;  // 确保 lc > 0
        for (auto& term : f)
            term.second /= c;
        return {std::move(c), std::move(f)};
    }

    // §7.4 子集乘积 mod m (带对称约化)
    inline upolynomial_<ZZ> __subset_product_mod(
        const std::vector<upolynomial_<ZZ>>& factors,
        const std::vector<size_t>& subset,
        const ZZ& lc_f,
        const ZZ& m)
    {
        upolynomial_<ZZ> prod;
        prod.push_back(std::make_pair(umonomial(0), lc_f));
        for (size_t idx : subset)
        {
            prod = prod * factors[idx];
            prod.normalization();
            __upoly_mod_coeff(prod, m);
        }
        return __upoly_symmetric_mod(prod, m);
    }

    // 获取 ZZ 多项式的常数项
    inline ZZ __upoly_const_term(const upolynomial_<ZZ>& f)
    {
        if (f.empty()) return ZZ(0);
        if (f.back().first.deg() == 0) return f.back().second;
        return ZZ(0);
    }

    // §7.5 Zassenhaus 因子重组
    inline std::vector<upolynomial_<ZZ>>
    __factor_recombine(
        const upolynomial_<ZZ>& f,
        const std::vector<upolynomial_<ZZ>>& lifted,
        const ZZ& m)
    {
        if (lifted.size() <= 1)
        {
            if (get_deg(f) > 0)
                return {f};
            return {};
        }

        size_t r = lifted.size();

        // T: 可用因子索引集合
        std::vector<int> T(r);
        std::iota(T.begin(), T.end(), 0);  // {0, 1, ..., r-1}
        upolynomial_<ZZ> f_star = f;
        std::vector<upolynomial_<ZZ>> result;

        // 字典序下一组合: idx[0..s-1] 是 [0, n) 中的 s 个下标
        // 返回 false 表示已穷举
        auto next_combination = [](std::vector<int>& idx, int n) -> bool {
            int s = (int)idx.size();
            int i = s - 1;
            while (i >= 0 && idx[i] == n - s + i)
                --i;
            if (i < 0) return false;
            ++idx[i];
            for (int j = i + 1; j < s; ++j)
                idx[j] = idx[j - 1] + 1;
            return true;
        };

        bool found;
        int s = 1;
        while (2 * s <= (int)T.size())
        {
            found = false;
            int n_active = (int)T.size();
            if (s > n_active / 2) break;

            // idx[0..s-1]: 当前子集在 T 中的下标
            std::vector<int> idx(s);
            std::iota(idx.begin(), idx.end(), 0);  // {0, 1, ..., s-1}

            do {
                // 映射回原始因子索引
                std::vector<size_t> S_idx(s);
                for (int j = 0; j < s; ++j)
                    S_idx[j] = T[idx[j]];

                ZZ lc_fstar = f_star.front().second;

                // Hensel 提升将 lc(f) 分配给了 lifted[0]。
                // 当子集包含 index 0 时，子集乘积已自带 lc，
                // 无需额外乘 lc_fstar。
                bool subset_has_lc = false;
                for (size_t i : S_idx)
                    if (i == 0) { subset_has_lc = true; break; }
                ZZ lc_mult = subset_has_lc ? ZZ(1) : lc_fstar;

                // === 剪枝 1: 首项系数检查 ===
                ZZ lc_prod = lc_mult;
                for (size_t i : S_idx)
                    lc_prod *= lifted[i].front().second;
                lc_prod = __symmetric_mod(lc_prod, m);
                ZZ lc_sq = lc_fstar * lc_fstar;
                if (lc_prod != ZZ(0))
                {
                    ZZ rem;
                    ZZ::fdiv_r(rem, lc_sq, lc_prod);
                    if (rem != ZZ(0))
                        continue;
                }

                {
                    // === 剪枝 2: 常数项检查 ===
                    ZZ c_prod = lc_mult;
                    for (size_t i : S_idx)
                        c_prod *= __upoly_const_term(lifted[i]);
                    c_prod = __symmetric_mod(c_prod, m);
                    ZZ fstar_const = lc_fstar * __upoly_const_term(f_star);
                    if (c_prod != ZZ(0))
                    {
                        ZZ rem2;
                        ZZ::fdiv_r(rem2, fstar_const, c_prod);
                        if (rem2 != ZZ(0))
                            continue;
                    }

                    // === 完整验证: 试除 ===
                    auto g = __subset_product_mod(lifted, S_idx, lc_mult, m);
                    auto [c_g, pp_g] = __upoly_primitive(std::move(g));

                    upolynomial_<ZZ> q_trial, r_trial;
                    pair_vec_div(q_trial.data(), r_trial.data(),
                                 f_star.data(), pp_g.data(), f_star.comp());

                    if (r_trial.empty())
                    {
                        // 找到真因子
                        result.push_back(std::move(pp_g));
                        auto [c_q, pp_q] = __upoly_primitive(std::move(q_trial));
                        f_star = std::move(pp_q);

                        // 从 T 中删除已用因子 (逆序删除保持下标有效)
                        for (int j = s - 1; j >= 0; --j)
                            T.erase(T.begin() + idx[j]);
                        found = true;
                        break;
                    }
                }
            } while (next_combination(idx, n_active));

            if (found)
            {
                s = 1;  // 重新从 s=1 开始
                continue;
            }
            ++s;
        }

        // 剩余的 f* 本身不可约
        if (!f_star.empty() && get_deg(f_star) > 0)
            result.push_back(std::move(f_star));

        // 按 degree 排序
        std::sort(result.begin(), result.end(),
            [](const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b) {
                return get_deg(a) < get_deg(b);
            });

        return result;
    }

    // ================================================================
    // §8.2 upolynomial → polynomial 转换辅助（委托 poly_convert）
    // ================================================================

    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> __upoly_to_poly(
        const upolynomial_<ZZ>& up,
        const variable& var,
        const lex_<var_order>* comp_ptr)
    {
        polynomial_<ZZ,lex_<var_order>> result(comp_ptr);
        poly_convert(up, result, var);
        return result;
    }

    // ================================================================
    // §8.3 素数选择
    // ================================================================

    struct __prime_selection_result {
        uint32_t prime;
        std::vector<upolynomial_<Zp>> factors;
        bool irreducible;
    };

    inline __prime_selection_result __select_prime(const upolynomial_<ZZ>& f)
    {
        assert(!f.empty() && get_deg(f) >= 2);
        ZZ lc_f = f.front().second;

        __prime_selection_result best;
        best.prime = 0;
        best.irreducible = false;
        size_t best_count = SIZE_MAX;

        int64_t deg_f = get_deg(f);
        size_t max_tries = 30;
        constexpr size_t PRIME_TABLE_SIZE = 9999;  // boost::math::prime 表上限
        std::mt19937 rng(42);

        for (size_t idx = 0, tried = 0; tried < max_tries; ++idx)
        {
            if (idx >= PRIME_TABLE_SIZE)
                throw std::runtime_error(
                    "factorize: exhausted prime table without finding a suitable prime");
            uint32_t p = boost::math::prime((unsigned)idx);

            // 跳过 lc(f) mod p == 0 的素数
            ZZ lc_mod;
            ZZ::fdiv_r(lc_mod, lc_f, ZZ(p));
            if (!lc_mod) continue;

            // f mod p
            auto fp = polynomial_mod(f, p);
            if (fp.empty() || get_deg(fp) != deg_f)
                continue;

            // 检查 f mod p 是否无平方: gcd(fp, fp') == 1
            auto fp_deriv = derivative(fp);
            if (fp_deriv.empty()) continue;
            auto g = polynomial_GCD(fp, fp_deriv);
            if (get_deg(g) > 0) continue;

            ++tried;

            // 已知 fp 无平方，直接 DDF + EDF 跳过 __squarefree_Zp
            __upoly_make_monic(fp);
            auto ddf = __ddf_Zp(fp);

            std::vector<upolynomial_<Zp>> irr_factors;
            for (auto& [gk, dk] : ddf)
            {
                std::vector<upolynomial_<Zp>> edf_out;
                __edf_Zp(edf_out, gk, dk, rng);
                for (auto& hi : edf_out)
                    irr_factors.push_back(std::move(hi));
            }

            size_t nfactors = irr_factors.size();

            // 只有 1 个因子 → 不可约
            if (nfactors <= 1)
            {
                best.prime = p;
                best.factors = std::move(irr_factors);
                if (best.factors.empty())
                    best.factors.push_back(std::move(fp));
                best.irreducible = true;
                return best;
            }

            if (nfactors < best_count)
            {
                best_count = nfactors;
                best.prime = p;
                best.factors = std::move(irr_factors);
            }

            // 因子数 > deg/2 时扩展尝试到 50 个
            if (best_count > (size_t)(deg_f / 2) && tried == 30)
                max_tries = 50;
        }

        best.irreducible = false;
        return best;
    }

    // ================================================================
    // §8.4 无平方 ZZ 因式分解
    // ================================================================

    // 前置条件: f 是无平方、本原、deg >= 2 的 ZZ[x] 多项式
    inline std::vector<upolynomial_<ZZ>>
    __factor_squarefree_primitive_ZZ(const upolynomial_<ZZ>& f)
    {
        assert(!f.empty() && get_deg(f) >= 2);

        // 选择素数
        auto sel = __select_prime(f);

        if (sel.irreducible || sel.factors.size() <= 1)
            return {f};

        // Hensel 提升
        auto [lifted, modulus] = __hensel_lift(f, sel.factors, sel.prime);

        // 因子重组
        return __factor_recombine(f, lifted, modulus);
    }

    // ================================================================
    // §8.8 factorize: upolynomial_ZZ
    // ================================================================

    inline factorization<upolynomial_<ZZ>>
    factorize(const upolynomial_<ZZ>& F)
    {
        factorization<upolynomial_<ZZ>> result;
        result.content = ZZ(1);

        if (F.empty())
        {
            result.content = ZZ(0);
            return result;
        }

        if (is_number(F))
        {
            result.content = F.front().second;
            return result;
        }

        // 提取内容，本原化
        auto [ct, uf_prim] = __upoly_primitive(upolynomial_<ZZ>(F));
        result.content = ct;

        if (get_deg(uf_prim) <= 1)
        {
            result.factors.push_back({std::move(uf_prim), 1});
            return result;
        }

        // 转 polynomial_ZZ 做无平方分解（squarefreefactorize 仅定义在 polynomial 上）
        variable __x("x");
        polynomial_<ZZ,lex> poly_prim;
        poly_convert(uf_prim, poly_prim, __x);
        auto sqf = squarefreefactorize(poly_prim);

        for (auto& [sqf_factor, mult] : sqf)
        {
            if (is_number(sqf_factor))
            {
                result.content *= sqf_factor.front().second;
                continue;
            }

            upolynomial_<ZZ> usqf;
            poly_convert(sqf_factor, usqf);

            if (get_deg(usqf) <= 1)
            {
                result.factors.push_back({std::move(usqf), mult});
                continue;
            }

            auto irr_factors = __factor_squarefree_primitive_ZZ(usqf);
            for (auto& irr : irr_factors)
            {
                if (irr.front().second < ZZ(0))
                {
                    for (auto& term : irr)
                        term.second = -term.second;
                }
                result.factors.push_back({std::move(irr), mult});
            }
        }

        std::sort(result.factors.begin(), result.factors.end(),
            [](const auto& a, const auto& b) {
                return degree(a.first) < degree(b.first);
            });

#ifndef NDEBUG
        {
            // 转换为 polynomial_<ZZ,lex> 验证乘积 == 原多项式
            variable __xdbg("x");
            polynomial_<ZZ,lex> check_poly;
            polynomial_<ZZ,lex> f_poly;
            poly_convert(F, f_poly, __xdbg);

            check_poly.data().push_back({{}, result.content});
            check_poly.normalization();
            for (const auto& [fi, ei] : result.factors)
            {
                polynomial_<ZZ,lex> fi_poly;
                poly_convert(fi, fi_poly, __xdbg);
                for (uint64_t e = 0; e < ei; ++e)
                {
                    check_poly = check_poly * fi_poly;
                    check_poly.normalization();
                }
            }
            assert(check_poly == f_poly);
        }
#endif

        return result;
    }

} // namespace clpoly

#endif // CLPOLY_POLYNOMIAL_FACTORIZE_UNIVAR_HH
