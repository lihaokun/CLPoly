/**
 * @file polynomial_factorize_zp.hh
 * @brief Zp[x] 单变量因式分解辅助函数与算法 (M1: squarefree / DDF / EDF)
 *
 * 本文件可独立包含。供 polynomial_factorize_univar.hh 使用。
 */
#ifndef CLPOLY_POLYNOMIAL_FACTORIZE_ZP_HH
#define CLPOLY_POLYNOMIAL_FACTORIZE_ZP_HH

#include <clpoly/upolynomial.hh>
#include <clpoly/polynomial_gcd.hh>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>

namespace clpoly{

    // Helper: 构造 Zp(int64_t, uint64_t) 避免重载歧义
    inline Zp __make_zp(int64_t val, uint64_t p) { return Zp(val, p); }

    // ================================================================
    // §4. 辅助函数层
    // ================================================================

    // §4.1.1 首一化
    inline Zp __upoly_make_monic(upolynomial_<Zp>& f)
    {
        assert(!f.empty());
        Zp lc = f.front().second;
        if (lc.number() == 1) return lc;
        Zp lc_inv = lc.inv();
        for (auto& term : f)
            term.second *= lc_inv;
        return lc;
    }

    // §4.1.2 多项式取模: f mod g in Zp[x]
    inline upolynomial_<Zp> __upoly_mod(
        const upolynomial_<Zp>& f,
        const upolynomial_<Zp>& g)
    {
        upolynomial_<Zp> q, r;
        pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
        return r;
    }

    // §4.1.3 商和余式: f = q*g + r in Zp[x]
    inline void __upoly_divmod(
        upolynomial_<Zp>& q,
        upolynomial_<Zp>& r,
        const upolynomial_<Zp>& f,
        const upolynomial_<Zp>& g)
    {
        pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
    }

    // §4.1.7 模幂: base^exp mod modpoly in Zp[x]
    inline upolynomial_<Zp> __upoly_powmod(
        const upolynomial_<Zp>& base,
        const ZZ& exp,
        const upolynomial_<Zp>& modpoly)
    {
        assert(!modpoly.empty());
        uint64_t p = modpoly.front().second.prime();
        Zp one = __make_zp(1, p);
        upolynomial_<Zp> result({std::make_pair(umonomial(0), one)});
        upolynomial_<Zp> b = __upoly_mod(base, modpoly);
        ZZ e = exp;
        while (e > 0)
        {
            if (e % 2 != 0)
            {
                result = result * b;
                result = __upoly_mod(result, modpoly);
            }
            e = e / 2;
            if (e > 0)
            {
                b = b * b;
                b = __upoly_mod(b, modpoly);
            }
        }
        return result;
    }

    // §4.1.8 随机多项式: 度数 < max_deg 的随机 Zp[x] 多项式
    inline upolynomial_<Zp> __upoly_random(
        int64_t max_deg,
        uint64_t p,
        std::mt19937& rng)
    {
        std::uniform_int_distribution<uint64_t> dist(0, (uint64_t)(p - 1));
        upolynomial_<Zp> result;
        for (int64_t d = max_deg - 1; d >= 0; --d)
        {
            uint64_t c = dist(rng);
            if (c != 0)
                result.push_back(std::make_pair(umonomial(d), Zp(c, p)));
        }
        return result;
    }

    // ================================================================
    // §5. M1: Zp 上的单变量分解
    // ================================================================

    // §5.3.1 p 次根提取: f(x) = g(x^p) => g
    inline upolynomial_<Zp> __extract_pth_root(const upolynomial_<Zp>& f)
    {
        uint64_t p = f.front().second.prime();
        upolynomial_<Zp> g;
        g.reserve(f.size());
        for (auto& term : f)
        {
            assert(term.first.deg() % p == 0);
            g.push_back(std::make_pair(umonomial(term.first.deg() / p), term.second));
        }
        return g;
    }

    // §5.3 Zp 上的无平方分解
    inline std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
    __squarefree_Zp(const upolynomial_<Zp>& f)
    {
        assert(!f.empty());
        std::vector<std::pair<upolynomial_<Zp>, uint64_t>> result;
        uint64_t p = f.front().second.prime();

        auto f_deriv = derivative(f);

        if (f_deriv.empty())
        {
            auto g = __extract_pth_root(f);
            __upoly_make_monic(g);
            auto sub = __squarefree_Zp(g);
            for (auto& si_ei : sub)
                result.push_back({std::move(si_ei.first), si_ei.second * p});
            return result;
        }

        auto c = polynomial_GCD(f, f_deriv);
        // w = f / c
        upolynomial_<Zp> w;
        pair_vec_div(w.data(), f.data(), c.data(), f.comp());
        w.normalization();

        uint64_t i = 1;
        while (!w.empty() && get_deg(w) > 0)
        {
            auto y = polynomial_GCD(w, c);
            // z = w / y
            upolynomial_<Zp> z;
            pair_vec_div(z.data(), w.data(), y.data(), w.comp());
            z.normalization();
            if (!z.empty() && get_deg(z) > 0)
            {
                __upoly_make_monic(z);
                result.push_back({std::move(z), i});
            }
            // c = c / y
            upolynomial_<Zp> c_new;
            pair_vec_div(c_new.data(), c.data(), y.data(), c.comp());
            c = std::move(c_new);
            c.normalization();
            w = std::move(y);
            ++i;
        }

        if (!c.empty() && get_deg(c) > 0)
        {
            auto g = __extract_pth_root(c);
            __upoly_make_monic(g);
            auto sub = __squarefree_Zp(g);
            for (auto& sj_ej : sub)
                result.push_back({std::move(sj_ej.first), sj_ej.second * p});
        }

        return result;
    }

    // 减去 x 的辅助函数: 从多项式中减去 monomial x
    inline upolynomial_<Zp> __upoly_subtract_x(const upolynomial_<Zp>& h, uint64_t p)
    {
        upolynomial_<Zp> result;
        result.reserve(h.size() + 1);
        bool inserted = false;
        for (auto& term : h)
        {
            if (!inserted && term.first.deg() < 1)
            {
                // 插入 -x 项 (即 (p-1)*x)
                Zp neg_one = __make_zp((int64_t)(p - 1), p);
                result.push_back(std::make_pair(umonomial(1), neg_one));
                inserted = true;
            }
            if (term.first.deg() == 1)
            {
                Zp new_c = term.second - __make_zp(1, p);
                if (new_c.number() != 0)
                    result.push_back(std::make_pair(umonomial(1), new_c));
                inserted = true;
                continue;
            }
            result.push_back(term);
        }
        if (!inserted)
        {
            Zp neg_one = __make_zp((int64_t)(p - 1), p);
            result.push_back(std::make_pair(umonomial(1), neg_one));
        }
        result.normalization();
        return result;
    }

    // 减去常数 1 的辅助函数
    inline upolynomial_<Zp> __upoly_subtract_one(const upolynomial_<Zp>& h, uint64_t p)
    {
        upolynomial_<Zp> result;
        result.reserve(h.size() + 1);
        bool found = false;
        for (auto& term : h)
        {
            if (term.first.deg() == 0)
            {
                Zp new_c = term.second - __make_zp(1, p);
                if (new_c.number() != 0)
                    result.push_back(std::make_pair(umonomial(0), new_c));
                found = true;
            }
            else
            {
                result.push_back(term);
            }
        }
        if (!found)
        {
            // 加入 -1 = p-1
            result.push_back(std::make_pair(umonomial(0), __make_zp((int64_t)(p - 1), p)));
            result.normalization();
        }
        return result;
    }

    // §5.4 DDF: 按度数分组
    inline std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
    __ddf_Zp(const upolynomial_<Zp>& f)
    {
        assert(!f.empty());
        uint64_t p = f.front().second.prime();
        std::vector<std::pair<upolynomial_<Zp>, uint64_t>> result;

        // h = x
        upolynomial_<Zp> h({std::make_pair(umonomial(1), __make_zp(1, p))});
        upolynomial_<Zp> f_star = f;

        for (uint64_t d = 1; ; ++d)
        {
            if (get_deg(f_star) < (int64_t)(2 * d))
                break;

            // h = h^p mod f*
            h = __upoly_powmod(h, ZZ(p), f_star);

            // h - x
            auto h_minus_x = __upoly_subtract_x(h, p);

            auto gd = polynomial_GCD(h_minus_x, f_star);

            if (!gd.empty() && get_deg(gd) > 0)
            {
                __upoly_make_monic(gd);
                result.push_back({gd, d});
                // f* = f* / gd
                upolynomial_<Zp> f_new;
                pair_vec_div(f_new.data(), f_star.data(), gd.data(), f_star.comp());
                f_star = std::move(f_new);
                f_star.normalization();
                // h = h mod f*
                h = __upoly_mod(h, f_star);
            }
        }

        if (!f_star.empty() && get_deg(f_star) > 0)
        {
            __upoly_make_monic(f_star);
            result.push_back({std::move(f_star), (uint64_t)get_deg(f_star)});
        }

        return result;
    }

    // §5.5 EDF: 等度分裂 (Cantor-Zassenhaus)
    inline void __edf_Zp(
        std::vector<upolynomial_<Zp>>& result,
        const upolynomial_<Zp>& f,
        uint64_t d,
        std::mt19937& rng)
    {
        if ((uint64_t)get_deg(f) == d)
        {
            auto f_copy = f;
            __upoly_make_monic(f_copy);
            result.push_back(std::move(f_copy));
            return;
        }
        if (get_deg(f) <= 0)
            return;

        uint64_t p = f.front().second.prime();
        int64_t n = get_deg(f);

        while (true)
        {
            auto r = __upoly_random(n, p, rng);
            if (r.empty())
                continue;

            upolynomial_<Zp> g;
            if (p == 2)
            {
                // 特征 2: trace map
                g = __upoly_mod(r, f);
                for (uint64_t i = 1; i < d; ++i)
                {
                    g = g * g + r;
                    g = __upoly_mod(g, f);
                }
                g = polynomial_GCD(g, f);
            }
            else
            {
                // 奇特征: g = gcd(r^{(p^d-1)/2} - 1, f)
                ZZ exp = (pow(ZZ(p), d) - 1) / 2;
                auto g_pow = __upoly_powmod(r, exp, f);
                auto g_pow_minus_1 = __upoly_subtract_one(g_pow, p);
                g = polynomial_GCD(g_pow_minus_1, f);
            }

            if (get_deg(g) > 0 && get_deg(g) < get_deg(f))
            {
                // 成功分裂
                upolynomial_<Zp> h_part;
                pair_vec_div(h_part.data(), f.data(), g.data(), f.comp());
                h_part.normalization();
                __upoly_make_monic(g);
                __upoly_make_monic(h_part);
                __edf_Zp(result, g, d, rng);
                __edf_Zp(result, h_part, d, rng);
                return;
            }
        }
    }

    // §5.6 Zp 上完整不可约分解
    inline std::pair<Zp, std::vector<std::pair<upolynomial_<Zp>, uint64_t>>>
    __factor_Zp(upolynomial_<Zp> f)
    {
        assert(!f.empty());
        uint64_t p = f.front().second.prime();

        if (get_deg(f) <= 0)
            return {f.front().second, {}};

        Zp lc = __upoly_make_monic(f);

        auto sqf = __squarefree_Zp(f);

        std::vector<std::pair<upolynomial_<Zp>, uint64_t>> result;
        std::mt19937 rng(42);

        for (auto& sj_ej : sqf)
        {
            auto ddf = __ddf_Zp(sj_ej.first);
            for (auto& gk_dk : ddf)
            {
                std::vector<upolynomial_<Zp>> factors_k;
                __edf_Zp(factors_k, gk_dk.first, gk_dk.second, rng);
                for (auto& hi : factors_k)
                    result.push_back({std::move(hi), sj_ej.second});
            }
        }

        std::sort(result.begin(), result.end(),
            [](const auto& a, const auto& b) {
                return get_deg(a.first) < get_deg(b.first);
            });

        return {lc, std::move(result)};
    }


} // namespace clpoly

#endif // CLPOLY_POLYNOMIAL_FACTORIZE_ZP_HH
