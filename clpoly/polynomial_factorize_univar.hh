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
#ifdef CLPOLY_PROFILE
#include <chrono>
#include <cstdio>
#include <atomic>
#endif

namespace clpoly{

#ifdef CLPOLY_PROFILE
    // ================================================================
    // Profiling accumulators (enabled only when -DCLPOLY_PROFILE)
    // ================================================================
    struct __profile_stats {
        std::atomic<long long> select_prime_ns{0};
        std::atomic<long long> hensel_lift_ns{0};
        std::atomic<long long> recombine_ns{0};
        std::atomic<long long> phase2_hensel_ns{0};
        std::atomic<long long> phase2_recombine_ns{0};
        std::atomic<long long> calls{0};
        std::atomic<long long> phase2_triggers{0};
    };
    inline __profile_stats __g_profile;

    inline void print_factorize_profile() {
        long long calls     = __g_profile.calls.load();
        long long sel_ns    = __g_profile.select_prime_ns.load();
        long long hen_ns    = __g_profile.hensel_lift_ns.load();
        long long rec_ns    = __g_profile.recombine_ns.load();
        long long h2_ns     = __g_profile.phase2_hensel_ns.load();
        long long r2_ns     = __g_profile.phase2_recombine_ns.load();
        long long total_ns  = sel_ns + hen_ns + rec_ns + h2_ns + r2_ns;
        long long ph2       = __g_profile.phase2_triggers.load();
        if (calls == 0 || total_ns == 0) { printf("[profile] no data\n"); return; }
        auto pct = [&](long long x) { return 100.0 * x / total_ns; };
        auto us  = [](long long ns) { return ns / 1000.0; };
        printf("\n===== factorize profile (%lld calls, %lld phase-2 triggers) =====\n", calls, ph2);
        printf("  select_prime       %9.1f us  %5.1f%%\n", us(sel_ns), pct(sel_ns));
        printf("  hensel_lift P1     %9.1f us  %5.1f%%\n", us(hen_ns), pct(hen_ns));
        printf("  recombine   P1     %9.1f us  %5.1f%%\n", us(rec_ns), pct(rec_ns));
        printf("  hensel_lift P2     %9.1f us  %5.1f%%\n", us(h2_ns),  pct(h2_ns));
        printf("  recombine   P2     %9.1f us  %5.1f%%\n", us(r2_ns),  pct(r2_ns));
        printf("  total              %9.1f us\n",           us(total_ns));
        printf("  avg per call       %9.1f us\n",           us(total_ns / calls));
        printf("================================================================\n\n");
    }
    inline void reset_factorize_profile() {
        __g_profile.select_prime_ns = 0;
        __g_profile.hensel_lift_ns  = 0;
        __g_profile.recombine_ns    = 0;
        __g_profile.phase2_hensel_ns= 0;
        __g_profile.phase2_recombine_ns = 0;
        __g_profile.calls           = 0;
        __g_profile.phase2_triggers = 0;
    }
    #define _PROF_NOW() std::chrono::steady_clock::now()
    #define _PROF_NS(a,b) std::chrono::duration_cast<std::chrono::nanoseconds>((b)-(a)).count()
#endif // CLPOLY_PROFILE

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
        uint64_t p,
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
        poly_convert(g_zp, nodes[parent_idx].g);
        poly_convert(h_zp, nodes[parent_idx].h);
        poly_convert(s_zp, nodes[parent_idx].s);
        poly_convert(t_zp, nodes[parent_idx].t);
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
        uint64_t p)
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
    // a_target == 0：提升到完整 Mignotte 精度（默认）
    // a_target >  0：提升到 p^a_target（停止条件 m > p^a_target - 1）
    inline std::pair<std::vector<upolynomial_<ZZ>>, ZZ>
    __hensel_lift(
        const upolynomial_<ZZ>& f,
        const std::vector<upolynomial_<Zp>>& factors,
        uint64_t p,
        int a_target = 0)
    {
        assert(factors.size() >= 2);
        assert(!f.empty());

        // 1. 确定提升精度
        ZZ target;
        if (a_target == 0) {
            ZZ B = __mignotte_bound(f);
            ZZ lc_f = f.front().second;
            if (lc_f < ZZ(0)) lc_f = -lc_f;
            target = ZZ(2) * lc_f * B;  // 需要 p^k > target
        } else {
            // target = p^a_target - 1：使 while(m <= target) 等价于 while(m < p^a_target)
            target = ZZ(1);
            for (int i = 0; i < a_target; ++i)
                target *= ZZ(p);
            target -= ZZ(1);
        }

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
    // §6.7 线性 Hensel 提升辅助 (P1b)
    // ================================================================

    // M4: 启发式起始精度（FLINT 公式 + Mignotte 上限）
    inline int __heuristic_starting_precision(
        const upolynomial_<ZZ>& f,
        int                     r,
        uint64_t                p)
    {
        // 1. FLINT 启发式 a_h（double 精度）
        double logp  = std::log((double)p);
        int    min_b = (int)ZZ(p).sizeinbase(2);
        int    N     = (int)f.size() - 1;
        double a_h_d = std::ceil(
            (2.5 * r + min_b) * std::log(2.0) / logp
            + std::log((double)(N + 1)) / (2.0 * logp));
        int a_h = std::max(1, (int)a_h_d);

        // 2. Mignotte 精度 a_mig（循环逼近）
        ZZ B_mig = __mignotte_bound(f);
        ZZ lc_f  = f.front().second;
        if (lc_f < ZZ(0)) lc_f = -lc_f;
        ZZ target = ZZ(2) * lc_f * B_mig;
        int a_mig = 0;
        ZZ  pa(1);
        while (pa <= target) { pa *= ZZ(p); ++a_mig; }

        return std::min(a_mig, a_h);
    }

    // M1: 二叉树单步线性 Hensel 提升（m → m·p）
    // 前置条件：node.g * node.h ≡ f (mod m)；node.s, node.t 固定 mod p
    // 后置条件：node.g * node.h ≡ f (mod m·p)；s, t 不变
    inline void __hensel_step_linear(
        __hensel_node&          node,
        const upolynomial_<ZZ>& f,
        const ZZ&               m,   // 当前模数 p^a
        uint64_t                p)   // 素数
    {
        ZZ p_zz(p);
        ZZ mp = m * p_zz;    // 新模数 p^(a+1)

        // 1. e = (f - g*h) / m  mod p
        upolynomial_<ZZ> gh = node.g * node.h;
        upolynomial_<ZZ> e  = f - gh;
        {
            auto it = e.data().begin(), out = it;
            for (; it != e.data().end(); ++it)
            {
                ZZ::fdiv_q(it->second, it->second, m);    // 精确整除
                ZZ::fdiv_r(it->second, it->second, p_zz); // mod p
                if (it->second) { if (out != it) *out = std::move(*it); ++out; }
            }
            e.data().erase(out, e.data().end());
        }
        if (e.empty()) return;  // 已满足精度，无需修正

        // 2. se = s * e，divmod by h (mod p) → σ = remainder
        upolynomial_<ZZ> se = node.s * e;
        upolynomial_<ZZ> q_se, sigma;
        __upoly_divmod_mod(q_se, sigma, se, node.h, p_zz);

        // 3. τ = t*e + q_se*g  (mod p)（包含商项，保证任意度 e 正确）
        upolynomial_<ZZ> tau = node.t * e + q_se * node.g;
        __upoly_mod_coeff(tau, p_zz);

        // 4. g' = g + m*τ (mod mp)
        for (auto& term : tau.data()) term.second *= m;
        node.g = node.g + tau;
        __upoly_mod_coeff(node.g, mp);

        // 5. h' = h + m*σ (mod mp)
        for (auto& term : sigma.data()) term.second *= m;
        node.h = node.h + sigma;
        __upoly_mod_coeff(node.h, mp);
        // 注：s, t 不更新（线性提升固定 Bézout 于 mod p）
    }

    // 线性 Hensel 树提升（单步，m → m·p，递归从顶向下）
    inline void __hensel_lift_linear_recursive(
        std::vector<__hensel_node>& nodes,
        int                         idx,
        const upolynomial_<ZZ>&     f,   // 当前节点目标
        const ZZ&                   m,
        uint64_t                    p)
    {
        __hensel_step_linear(nodes[idx], f, m, p);
        if (nodes[idx].left >= 0)
            __hensel_lift_linear_recursive(nodes, nodes[idx].left, nodes[idx].g, m, p);
        if (nodes[idx].right >= 0)
            __hensel_lift_linear_recursive(nodes, nodes[idx].right, nodes[idx].h, m, p);
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
    __zassenhaus_recombine(
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
    // §7.6 van Hoeij LLL 重组：辅助常量与类型
    // ================================================================

    // r ≤ ZASSENHAUS_THRESHOLD 时使用 Zassenhaus，否则用 van Hoeij LLL
    static constexpr int ZASSENHAUS_THRESHOLD = 10;

    // 格矩阵：行主序，M[i] 是第 i 个格基向量（ZZ 整数坐标）
    using LLLMatrix = std::vector<std::vector<ZZ>>;

    // ================================================================
    // §7.7 M1: __cld_polys
    //   计算 CLD 多项式：C_i = (f_star / h_i) · h_i'  mod m，对称约化
    //   前置条件：每个 h_i 首一且精确整除 f_star (mod m)
    // ================================================================

    inline std::vector<upolynomial_<ZZ>>
    __cld_polys(
        const upolynomial_<ZZ>&              f_star,
        const std::vector<upolynomial_<ZZ>>& active_factors,
        const ZZ&                            m)
    {
        std::vector<upolynomial_<ZZ>> result;
        result.reserve(active_factors.size());

        for (const auto& h_i : active_factors)
        {
            // 步骤 1：精确除法 q_i = f_star / h_i  in Z_m[x]
            upolynomial_<ZZ> q_i, r_i;
            __upoly_divmod_mod(q_i, r_i, f_star, h_i, m);
            assert(r_i.empty());  // Hensel 提升保证 h_i | f_star (mod m)

            // 步骤 2：形式导数 h_i'，系数 mod m
            auto h_prime = derivative(h_i);
            __upoly_mod_coeff(h_prime, m);

            // 步骤 3：C_i = q_i · h_i'  mod m，再对称约化到 (-m/2, m/2]
            auto C_i = __upoly_mul_mod(q_i, h_prime, m);
            C_i      = __upoly_symmetric_mod(C_i, m);

            result.push_back(std::move(C_i));
        }
        return result;
    }

    // ================================================================
    // §7.8 M2: __build_cld_matrix
    //   按内螺旋顺序将 CLD 系数列喂入格矩阵 M，扩展 M 直到加入 J_target 新列。
    //   前置条件：M 已是 (r+J_cur)×(r+J_cur) 的整数方阵（行主序）
    //   后置条件：M 扩展为 (r+J_cur+J_new)×(r+J_cur+J_new)
    //   注：初版不做列过滤，所有螺旋位置均接受，后续可加过滤条件优化
    // ================================================================

    inline int __build_cld_matrix(
        LLLMatrix&                           M,
        const std::vector<upolynomial_<ZZ>>& cld,
        int                                  J_cur,
        int                                  J_target,
        const ZZ&                            /*m*/)
    {
        int r = (int)cld.size();

        // 辅助：按次数取 CLD 多项式的系数（upolynomial_ 不支持下标索引）
        auto upoly_coeff = [](const upolynomial_<ZZ>& p, int deg) -> ZZ {
            for (const auto& term : p)
                if ((int)term.first.deg() == deg) return term.second;
            return ZZ(0);
        };

        // N = 螺旋总位置数 = max(deg(C_i)) + 1（CLD 次数 ≤ deg(f)-1，故 N ≤ deg(f)）
        int N = 0;
        for (const auto& c : cld)
            if (!c.empty())
                N = std::max(N, (int)c.front().first.deg() + 1);

        int J_new = 0;
        for (int k = J_cur; J_new < J_target && k < N; ++k)
        {
            // 螺旋位置 k → CLD 系数次数 col_idx
            int col_idx = (k % 2 == 0) ? (k / 2) : (N - 1 - (k - 1) / 2);

            // j = 加入本列前 M 中已有的 J 列数（J_cur + J_new）
            int j = J_cur + J_new;

            // 步骤 1：扩展现有行，各追加一个 0
            for (auto& row : M)
                row.push_back(ZZ(0));

            // 步骤 2：新增数据行，长度 r + j + 1
            //   前 r 项 = cld[i][col_idx]，中间 j 项 = 0，末项 = 1（对角 identity）
            std::vector<ZZ> new_row(r + j + 1, ZZ(0));
            for (int i = 0; i < r; ++i)
                new_row[i] = upoly_coeff(cld[i], col_idx);
            new_row[r + j] = ZZ(1);

            M.push_back(std::move(new_row));
            ++J_new;
        }
        return J_new;
    }

    // ================================================================
    // §7.9 M3: __lll_reduce
    //   整数 LLL 基规约（Cohen §2.6 算法 2.6.3，δ=3/4），记录幺模变换矩阵 U。
    //   前置条件：M 是 n×n 整数方阵（n ≥ 1）
    //   后置条件：M 已 LLL 规约；U 满足 M_new = U · M_old；返回 ‖M[i]‖² ≤ B 的行下标
    // ================================================================

    inline std::vector<int>
    __lll_reduce(
        LLLMatrix& M,
        LLLMatrix& U,
        const ZZ&  B)
    {
        int n = (int)M.size();

        // Gram-Schmidt 系数 mu[i][j] (i > j) 和平方范数 B_gs[i] = ‖b_i^*‖²
        std::vector<std::vector<QQ>> mu(n, std::vector<QQ>(n, QQ(0)));
        std::vector<QQ>              B_gs(n, QQ(0));

        // U 初始化为单位矩阵
        U.assign(n, std::vector<ZZ>(n, ZZ(0)));
        for (int i = 0; i < n; ++i)
            U[i][i] = ZZ(1);

        // 辅助：ZZ 精确内积（行向量）
        auto dot = [&](const std::vector<ZZ>& a, const std::vector<ZZ>& b) -> ZZ {
            ZZ s(0);
            for (int k = 0; k < (int)a.size(); ++k)
                s += a[k] * b[k];
            return s;
        };

        // 辅助：四舍五入到最近整数（ties 向上）：floor(q + 1/2)
        // = floor((2*num + den) / (2*den))（den > 0 保证）
        auto round_qq = [](const QQ& q) -> ZZ {
            ZZ a = q.get_num() * ZZ(2) + q.get_den();
            ZZ b = q.get_den() * ZZ(2);
            ZZ result;
            ZZ::fdiv_q(result, a, b);
            return result;
        };

        // 辅助：行操作（同步作用于 M 和 U）：M[i] -= c * M[j], U[i] -= c * U[j]
        auto row_sub = [&](int i, int j, const ZZ& c) {
            for (int k = 0; k < n; ++k) {
                M[i][k] -= c * M[j][k];
                U[i][k] -= c * U[j][k];
            }
        };
        auto row_swap = [&](int i, int j) {
            std::swap(M[i], M[j]);
            std::swap(U[i], U[j]);
        };

        // 初始化 Gram-Schmidt
        B_gs[0] = QQ(dot(M[0], M[0]), ZZ(1));
        for (int i = 1; i < n; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                // num = <b_i, b_j> - sum_{l<j} mu[i][l] * mu[j][l] * B_gs[l]
                //     = <b_i, b_j^*>
                QQ num(dot(M[i], M[j]), ZZ(1));
                for (int l = 0; l < j; ++l)
                    num -= mu[i][l] * mu[j][l] * B_gs[l];
                mu[i][j] = (B_gs[j] == QQ(0)) ? QQ(0) : num / B_gs[j];
            }
            // B_gs[i] = ‖b_i‖² - sum_{j<i} mu[i][j]² * B_gs[j] = ‖b_i^*‖²
            B_gs[i] = QQ(dot(M[i], M[i]), ZZ(1));
            for (int j = 0; j < i; ++j)
                B_gs[i] -= mu[i][j] * mu[i][j] * B_gs[j];
        }

        // LLL 主循环（Cohen §2.6 算法 2.6.3）
        int k = 1;
        while (k < n)
        {
            // 大小规约（j = k-1）
            ZZ q = round_qq(mu[k][k - 1]);
            if (q != ZZ(0))
            {
                row_sub(k, k - 1, q);
                for (int j = 0; j < k - 1; ++j)
                    mu[k][j] -= QQ(q, ZZ(1)) * mu[k - 1][j];
                mu[k][k - 1] -= QQ(q, ZZ(1));
            }

            // Lovász 条件：B_gs[k] >= (3/4 - mu[k][k-1]^2) * B_gs[k-1]
            QQ lhs = B_gs[k];
            QQ rhs = (QQ(3, 4) - mu[k][k - 1] * mu[k][k - 1]) * B_gs[k - 1];

            if (lhs >= rhs)
            {
                // 额外大小规约（j = k-2 down to 0）
                for (int j = k - 2; j >= 0; --j)
                {
                    ZZ q2 = round_qq(mu[k][j]);
                    if (q2 != ZZ(0))
                    {
                        row_sub(k, j, q2);
                        for (int l = 0; l < j; ++l)
                            mu[k][l] -= QQ(q2, ZZ(1)) * mu[j][l];
                        mu[k][j] -= QQ(q2, ZZ(1));
                    }
                }
                ++k;
            }
            else
            {
                // Lovász 条件失败：交换 M[k] 和 M[k-1]，更新 Gram-Schmidt
                QQ mu_old = mu[k][k - 1];
                QQ B_new  = B_gs[k] + mu_old * mu_old * B_gs[k - 1];
                QQ mu_new = (B_new != QQ(0)) ? mu_old * B_gs[k - 1] / B_new : QQ(0);
                if (B_new != QQ(0))
                {
                    B_gs[k]     = B_gs[k] * B_gs[k - 1] / B_new;
                    B_gs[k - 1] = B_new;
                }
                row_swap(k, k - 1);
                std::swap(mu[k], mu[k - 1]);
                // swap 后 mu[k][k-1] 为 0（来自原 mu[k-1][k-1]），显式写入
                mu[k][k - 1] = mu_new;
                // 修正 j > k 行的 mu[j][k] 和 mu[j][k-1]
                for (int j = k + 1; j < n; ++j)
                {
                    QQ t        = mu[j][k];
                    mu[j][k]    = mu[j][k - 1] - mu_old * t;
                    mu[j][k - 1] = t + mu_new * mu[j][k];  // mu[j][k] 已是更新后的值
                }
                k = std::max(k - 1, 1);
            }
        }

        // 收集满足 ‖M[i]‖² ≤ B 的行，按范数²升序返回
        std::vector<int> short_rows;
        for (int i = 0; i < n; ++i)
        {
            ZZ norm_sq = dot(M[i], M[i]);
            if (norm_sq <= B)
                short_rows.push_back(i);
        }
        std::sort(short_rows.begin(), short_rows.end(),
            [&](int a, int b) { return dot(M[a], M[a]) < dot(M[b], M[b]); });
        return short_rows;
    }

    // ================================================================
    // §7.10 M4: __extract_candidates
    //   从 LLL 幺模矩阵 U 的短向量行中提取候选因子子集（列等价类分组）。
    //   前置条件：short_rows 均为 U 的有效行下标；r ≤ n（U 为 n×n）
    //   后置条件：返回 active 因子下标的等价类分组（每类下标集合）
    // ================================================================

    inline std::vector<std::vector<int>>
    __extract_candidates(
        const std::vector<int>& short_rows,
        const LLLMatrix&        U,
        int                     r)
    {
        if (short_rows.empty()) return {};

        int s = (int)short_rows.size();

        // U_short[k][j] = U[short_rows[k]][j]，仅取前 r 列
        std::vector<std::vector<ZZ>> U_short(s, std::vector<ZZ>(r));
        for (int k = 0; k < s; ++k)
            for (int j = 0; j < r; ++j)
                U_short[k][j] = U[short_rows[k]][j];

        // 列等价类分组：part[j] = j 所在等价类编号
        // 两列相等 ⟺ 对所有 k：U_short[k][j1] == U_short[k][j2]
        std::vector<int> part(r, -1);
        int num_classes = 0;
        for (int j = 0; j < r; ++j)
        {
            if (part[j] != -1) continue;
            part[j] = num_classes;
            for (int j2 = j + 1; j2 < r; ++j2)
            {
                if (part[j2] != -1) continue;
                bool equal = true;
                for (int k = 0; k < s; ++k)
                    if (U_short[k][j] != U_short[k][j2]) { equal = false; break; }
                if (equal)
                    part[j2] = num_classes;
            }
            ++num_classes;
        }

        // 按等价类编号收集候选子集
        std::vector<std::vector<int>> candidates(num_classes);
        for (int j = 0; j < r; ++j)
            candidates[part[j]].push_back(j);

        return candidates;
    }

    // ================================================================
    // §7.11 M5: __vanhoeij_recombine
    //   van Hoeij LLL 因子重组主控循环。
    //   前置条件：f 本原，deg(f)≥2，lc(f)>0；|lifted|≥2；m > 2·lc(f)·B_Mig(f)
    // ================================================================

    inline std::vector<upolynomial_<ZZ>>
    __vanhoeij_recombine(
        const upolynomial_<ZZ>&              f,
        const std::vector<upolynomial_<ZZ>>& lifted,
        const ZZ&                            m)
    {
        int r = (int)lifted.size();
        int N = (int)get_deg(f);

        // === 初始化参数 ===
        int U_exp = (int)ZZ(r > 20 ? r : 20).sizeinbase(2);
        ZZ  B     = ZZ(r + 1) * (ZZ(1) << (2 * U_exp));
        int J_max = (N + 1) / 2;
        int J0    = (3 * r > N + 1) ? 30 : 10;
        J0        = std::min(J0, J_max);

        // 活跃因子下标集合：active[k] = lifted 中的原始下标
        std::vector<int> active(r);
        std::iota(active.begin(), active.end(), 0);

        upolynomial_<ZZ>              f_star = f;
        std::vector<upolynomial_<ZZ>> result;

        // 构造初始格矩阵辅助 lambda
        auto make_initial_M = [](int rr, int U_exp_) -> LLLMatrix {
            LLLMatrix M(rr, std::vector<ZZ>(rr, ZZ(0)));
            ZZ scale = ZZ(1) << U_exp_;
            for (int i = 0; i < rr; ++i)
                M[i][i] = scale;
            return M;
        };

        LLLMatrix M = make_initial_M(r, U_exp);
        int J_cur    = 0;
        // 先以 J_target=0 做一次对角 LLL（等价于 s=1 Zassenhaus，近零开销）：
        // 若所有模因子已各自对应真因子，可在不建 CLD 列的情况下提取所有因子；
        // 若对角 LLL 未提取全部因子，再增加 CLD 列重试（J_target → J0 → 2·J0 → …）。
        int J_target = 0;

        while ((int)active.size() > 1)
        {
            // [M1] 计算当前活跃因子的 CLD 多项式（仅在需要添加列时才计算）
            std::vector<upolynomial_<ZZ>> active_lifted;
            for (int k : active)
                active_lifted.push_back(lifted[k]);
            std::vector<upolynomial_<ZZ>> cld;
            if (J_target > 0)
                cld = __cld_polys(f_star, active_lifted, m);

            // [M2] 喂入 CLD 列（J_target=0 时跳过，矩阵保持纯对角）
            int J_new = 0;
            if (J_target > 0)
            {
                J_new = __build_cld_matrix(M, cld, J_cur, J_target, m);
                J_cur += J_new;
            }

            // [M3] LLL 规约
            LLLMatrix U;
            auto short_rows = __lll_reduce(M, U, B);

            // [M4] 提取候选子集
            auto candidates = __extract_candidates(short_rows, U, (int)active.size());

            // === 候选子集验证（本轮内批量提取所有因子）===
            // consumed[j] = true 表示 active_lifted[j] 已在本轮被某候选消耗
            bool found_any = false;
            std::vector<bool> consumed((int)active.size(), false);
            int remaining_active = (int)active.size();

            for (auto& cand : candidates)
            {
                if (cand.empty()) continue;
                // 跳过覆盖全部剩余活跃因子的候选（等价于平凡分割）
                if ((int)cand.size() >= remaining_active) continue;
                // 跳过已被前一候选消耗的因子
                bool any_consumed = false;
                for (int k : cand)
                    if (consumed[k]) { any_consumed = true; break; }
                if (any_consumed) continue;

                // 构造候选乘积：lc(f_star) × ∏ active_lifted[k]，逐步 normalization + mod
                ZZ lc_fstar = f_star.front().second;
                upolynomial_<ZZ> g_trial;
                g_trial.push_back({umonomial(0), lc_fstar});
                for (int k : cand)
                {
                    g_trial = g_trial * active_lifted[k];
                    g_trial.normalization();
                    __upoly_mod_coeff(g_trial, m);
                }
                g_trial = __upoly_symmetric_mod(g_trial, m);

                auto [c_g, pp_g] = __upoly_primitive(std::move(g_trial));

                // 试除（针对当前 f_star，已从前一成功候选中约去）
                upolynomial_<ZZ> q_trial, r_trial;
                pair_vec_div(q_trial.data(), r_trial.data(),
                             f_star.data(), pp_g.data(), f_star.comp());

                if (!r_trial.empty()) continue;

                // 成功：记录因子，更新 f_star 和消耗标记
                result.push_back(std::move(pp_g));
                auto [c_q, pp_q] = __upoly_primitive(std::move(q_trial));
                f_star = std::move(pp_q);

                for (int k : cand) consumed[k] = true;
                remaining_active -= (int)cand.size();
                found_any = true;
            }

            if (found_any)
            {
                // 按逆序从 active 中移除本轮消耗的所有因子
                for (int j = (int)active.size() - 1; j >= 0; --j)
                    if (consumed[j]) active.erase(active.begin() + j);

                // 重建格矩阵参数
                int r_new   = (int)active.size();
                int U_exp_n = (int)ZZ(r_new > 20 ? r_new : 20).sizeinbase(2);
                B           = ZZ(r_new + 1) * (ZZ(1) << (2 * U_exp_n));
                M           = make_initial_M(r_new, U_exp_n);
                J_cur       = 0;
                J_target    = 0;  // 找到因子后重置：新子问题先尝试对角 LLL
                continue;
            }

            // 本轮无新因子：
            //   J_target == 0 时（对角 LLL 失败）：升至 J0，开始加入 CLD 列
            //   J_target > 0 时：翻倍 J_target，超出 J_max 则回退 Zassenhaus
            if (J_target == 0)
                J_target = J0;
            else
                J_target *= 2;
            if (J_target > J_max)
            {
                auto zass = __zassenhaus_recombine(f_star, active_lifted, m);
                for (auto& g : zass)
                    result.push_back(std::move(g));
                return result;
            }
        }

        // active 剩余 ≤ 1：若 f_star 仍有正次数则加入结果
        if (!f_star.empty() && get_deg(f_star) > 0)
            result.push_back(std::move(f_star));

        std::sort(result.begin(), result.end(),
            [](const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b) {
                return get_deg(a) < get_deg(b);
            });
        return result;
    }

    // ================================================================
    // §7.12 路由层 __factor_recombine
    //   lifted.size() ≤ ZASSENHAUS_THRESHOLD → Zassenhaus（O(2^r) 但 r 小）
    //   lifted.size() >  ZASSENHAUS_THRESHOLD → van Hoeij LLL（多项式复杂度）
    // ================================================================

    inline std::vector<upolynomial_<ZZ>>
    __factor_recombine(
        const upolynomial_<ZZ>&              f,
        const std::vector<upolynomial_<ZZ>>& lifted,
        const ZZ&                            m)
    {
        if ((int)lifted.size() <= ZASSENHAUS_THRESHOLD)
            return __zassenhaus_recombine(f, lifted, m);
        else
            return __vanhoeij_recombine(f, lifted, m);
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
        uint64_t prime;
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
        size_t max_tries = 3;
        constexpr size_t PRIME_TABLE_SIZE = 9999;  // boost::math::prime 表上限
        std::mt19937 rng(42);

        for (size_t idx = 0, tried = 0; tried < max_tries; ++idx)
        {
            if (idx >= PRIME_TABLE_SIZE)
                throw std::runtime_error(
                    "factorize: exhausted prime table without finding a suitable prime");
            uint64_t p = boost::math::prime((unsigned)idx);

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

        }

        best.irreducible = false;
        return best;
    }

    // ================================================================
    // §8.3 van Hoeij LLL 因子重组（P1a+P1b 统一入口）
    // 初始提升：二次 Hensel（快速到达完整精度）
    // 重组：van Hoeij LLL（所有 r，移除 Zassenhaus 阈值）
    // 安全网：Zassenhaus（当 LLL 列耗尽时）
    //
    // 注：线性 Hensel 基础设施（__hensel_step_linear、__hensel_lift_linear_recursive）
    //   已在 §6.7 实现，保留供未来真正的线性交织优化使用（P1b 完整版）。
    //   Van Hoeij LLL 需要在满精度（a_mig 或 a_h）下运行才能高效收敛；
    //   低精度下运行 LLL 会导致格基规约需 O(n^4) 次迭代而非 O(n)，性能急剧下降。
    // ================================================================

    inline std::vector<upolynomial_<ZZ>>
    __lll_factorize(
        const upolynomial_<ZZ>&              f,
        const std::vector<upolynomial_<Zp>>& factors,
        uint64_t                             p)
    {
        int r   = (int)factors.size();
        int a_h = __heuristic_starting_precision(f, r, p);

        // Phase 1：提升到启发式精度 a_h（当 a_h < a_mig 时节省提升代价）
#ifdef CLPOLY_PROFILE
        auto _t0 = _PROF_NOW();
#endif
        auto [lifted_h, m_h] = __hensel_lift(f, factors, p, a_h);
#ifdef CLPOLY_PROFILE
        auto _t1 = _PROF_NOW();
#endif
        auto result = __vanhoeij_recombine(f, lifted_h, m_h);
#ifdef CLPOLY_PROFILE
        auto _t2 = _PROF_NOW();
        __g_profile.hensel_lift_ns += _PROF_NS(_t0, _t1);
        __g_profile.recombine_ns   += _PROF_NS(_t1, _t2);
#endif

        // Phase 1 返回单个因子且 r > 1，说明精度不足（Mignotte 条件未满足）。
        // Phase 2：提升到完整 Mignotte 精度，保证正确性。
        if ((int)result.size() == 1 && r > 1) {
#ifdef CLPOLY_PROFILE
            auto _t3 = _PROF_NOW();
#endif
            auto [lifted_mig, m_mig] = __hensel_lift(f, factors, p);
#ifdef CLPOLY_PROFILE
            auto _t4 = _PROF_NOW();
#endif
            auto result2 = __vanhoeij_recombine(f, lifted_mig, m_mig);
#ifdef CLPOLY_PROFILE
            auto _t5 = _PROF_NOW();
            __g_profile.phase2_hensel_ns    += _PROF_NS(_t3, _t4);
            __g_profile.phase2_recombine_ns += _PROF_NS(_t4, _t5);
            __g_profile.phase2_triggers++;
#endif
            return result2;
        }
        return result;
    }

    // ================================================================
    // §8.4 无平方 ZZ 因式分解
    // ================================================================

    // 前置条件: f 是无平方、本原、deg >= 2 的 ZZ[x] 多项式
    inline std::vector<upolynomial_<ZZ>>
    __factor_squarefree_primitive_ZZ(const upolynomial_<ZZ>& f)
    {
        assert(!f.empty() && get_deg(f) >= 2);

#ifdef CLPOLY_PROFILE
        __g_profile.calls++;
        auto _ts0 = _PROF_NOW();
#endif
        // 选择素数
        auto sel = __select_prime(f);
#ifdef CLPOLY_PROFILE
        __g_profile.select_prime_ns += _PROF_NS(_ts0, _PROF_NOW());
#endif

        if (sel.irreducible || sel.factors.size() <= 1)
            return {f};

        // 二次 Hensel 提升 + van Hoeij LLL（移除 Zassenhaus 阈值，全 r 走 LLL）
        return __lll_factorize(f, sel.factors, sel.prime);
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
