/**
 * @file groebner.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief Gröbner basis computation
 */
#ifndef CLPOLY_GROEBNER_HH
#define CLPOLY_GROEBNER_HH

#include <clpoly/polynomial.hh>
#include <clpoly/polynomial_gcd.hh>
#include <clpoly/polynomial_convert.hh>
#include <vector>
#include <algorithm>
#include <cstdint>

namespace clpoly {

    // ===== M2: S-多项式 =====

    template<class Tc, class comp>
    inline polynomial_<Tc, comp> __s_polynomial(
        const polynomial_<Tc, comp>& f,
        const polynomial_<Tc, comp>& g)
    {
        assert(!f.empty() && !g.empty());
        assert(comp_consistent(f.comp(), g.comp()));

        const auto& LM_f = f.front().first;
        const auto& LC_f = f.front().second;
        const auto& LM_g = g.front().first;
        const auto& LC_g = g.front().second;

        auto L = lcm(LM_f, LM_g);
        auto m_f = L / LM_f;
        auto m_g = L / LM_g;

        // (1/LC_f) * m_f * f
        polynomial_<Tc, comp> poly_f(f);
        for (auto& term : poly_f.data())
        {
            term.first = term.first * m_f;
            term.second = term.second / LC_f;
        }

        // (1/LC_g) * m_g * g
        polynomial_<Tc, comp> poly_g(g);
        for (auto& term : poly_g.data())
        {
            term.first = term.first * m_g;
            term.second = term.second / LC_g;
        }

        // S = poly_f - poly_g
        polynomial_<Tc, comp> S(f.comp_ptr());
        pair_vec_sub(S.data(), poly_f.data(), poly_g.data(), f.comp());
        return S;
    }

    // 公开 API
    template<class Tc, class comp>
    inline polynomial_<Tc, comp> s_polynomial(
        const polynomial_<Tc, comp>& f,
        const polynomial_<Tc, comp>& g)
    {
        return __s_polynomial(f, g);
    }

    // ===== M3: Normal Form（多除数约化） =====

    template<class Tc, class comp>
    inline polynomial_<Tc, comp> __normal_form(
        const polynomial_<Tc, comp>& f,
        const std::vector<polynomial_<Tc, comp>>& G)
    {
        if (f.empty() || G.empty()) return f;

        polynomial_<Tc, comp> h(f);
        polynomial_<Tc, comp> r(f.comp_ptr());

        while (!h.empty())
        {
            bool divided = false;
            for (size_t i = 0; i < G.size(); ++i)
            {
                basic_monomial<comp> quot_m;
                if (is_divexact(quot_m, h.front().first, G[i].front().first))
                {
                    Tc c = h.front().second / G[i].front().second;
                    // h <- h - c * quot_m * G[i]
                    polynomial_<Tc, comp> tmp(G[i]);
                    for (auto& term : tmp.data())
                    {
                        term.first = term.first * quot_m;
                        term.second = term.second * c;
                    }
                    polynomial_<Tc, comp> h_new(h.comp_ptr());
                    pair_vec_sub(h_new.data(), h.data(), tmp.data(), h.comp());
                    h = std::move(h_new);
                    divided = true;
                    break;
                }
            }
            if (!divided)
            {
                r.data().push_back(h.front());
                h.data().erase(h.data().begin());
            }
        }
        return r;
    }

    template<class Tc, class comp>
    inline polynomial_<Tc, comp> __normal_form_with_sugar(
        const polynomial_<Tc, comp>& f,
        const std::vector<polynomial_<Tc, comp>>& G,
        const std::vector<int64_t>& sugar_vec,
        int64_t& f_sugar)
    {
        if (f.empty() || G.empty()) return f;

        polynomial_<Tc, comp> h(f);
        polynomial_<Tc, comp> r(f.comp_ptr());

        while (!h.empty())
        {
            bool divided = false;
            for (size_t i = 0; i < G.size(); ++i)
            {
                basic_monomial<comp> quot_m;
                if (is_divexact(quot_m, h.front().first, G[i].front().first))
                {
                    int64_t lm_h_deg = h.front().first.deg();
                    Tc c = h.front().second / G[i].front().second;
                    // h <- h - c * quot_m * G[i]
                    polynomial_<Tc, comp> tmp(G[i]);
                    for (auto& term : tmp.data())
                    {
                        term.first = term.first * quot_m;
                        term.second = term.second * c;
                    }
                    polynomial_<Tc, comp> h_new(h.comp_ptr());
                    pair_vec_sub(h_new.data(), h.data(), tmp.data(), h.comp());
                    h = std::move(h_new);
                    // sugar update
                    f_sugar = std::max(f_sugar,
                        sugar_vec[i] + lm_h_deg - G[i].front().first.deg());
                    divided = true;
                    break;
                }
            }
            if (!divided)
            {
                r.data().push_back(h.front());
                h.data().erase(h.data().begin());
            }
        }
        return r;
    }

    // 公开 API
    template<class Tc, class comp>
    inline polynomial_<Tc, comp> normal_form(
        const polynomial_<Tc, comp>& f,
        const std::vector<polynomial_<Tc, comp>>& G)
    {
        return __normal_form(f, G);
    }

    // ===== 辅助: 首一化（逐系数除以 LC） =====

    template<class Tc, class comp>
    inline void __make_monic(polynomial_<Tc, comp>& f)
    {
        if (f.empty()) return;
        Tc lc = f.front().second;
        if (lc == Tc(1)) return;
        for (auto& term : f.data())
            term.second = term.second / lc;
    }

    // ===== M4: 临界对管理 =====

    struct __critical_pair {
        size_t i, j;
        int64_t sugar;
        int64_t lcm_deg;
    };

    struct __cp_compare {
        bool operator()(const __critical_pair& a, const __critical_pair& b) const {
            if (a.sugar != b.sugar) return a.sugar < b.sugar;
            return a.lcm_deg < b.lcm_deg;
        }
    };

    template<class Tc, class comp>
    inline void __update(
        std::vector<polynomial_<Tc, comp>>& G,
        std::vector<__critical_pair>& pairs,
        const polynomial_<Tc, comp>& h,
        std::vector<int64_t>& sugar_vec,
        int64_t h_sugar)
    {
        size_t new_idx = G.size();
        const auto& LM_h = h.front().first;

        // 1. 生成新对候选
        std::vector<__critical_pair> new_pairs;
        for (size_t k = 0; k < G.size(); ++k)
        {
            auto L = lcm(G[k].front().first, LM_h);
            int64_t s = std::max(
                sugar_vec[k] + L.deg() - G[k].front().first.deg(),
                h_sugar + L.deg() - LM_h.deg());
            new_pairs.push_back({k, new_idx, s, L.deg()});
        }

        // 2. LCM 准则: 用 h 淘汰旧对
        std::vector<__critical_pair> filtered_pairs;
        for (auto& cp : pairs)
        {
            auto L_ij = lcm(G[cp.i].front().first, G[cp.j].front().first);
            basic_monomial<comp> quot_m;
            if (is_divexact(quot_m, L_ij, LM_h))
            {
                auto L_ih = lcm(G[cp.i].front().first, LM_h);
                auto L_hj = lcm(LM_h, G[cp.j].front().first);
                if (L_ih != L_ij && L_hj != L_ij)
                    continue; // 淘汰
            }
            filtered_pairs.push_back(cp);
        }
        pairs = std::move(filtered_pairs);

        // 3. 在新对中自身淘汰（保留 lcm 最小的）
        std::sort(new_pairs.begin(), new_pairs.end(),
            [](const __critical_pair& a, const __critical_pair& b) {
                return a.lcm_deg < b.lcm_deg;
            });

        // 预计算各新对的 lcm
        std::vector<basic_monomial<comp>> new_lcms;
        new_lcms.reserve(new_pairs.size());
        for (auto& p : new_pairs)
            new_lcms.push_back(lcm(G[p.i].front().first, LM_h));

        std::vector<__critical_pair> minimal_new;
        std::vector<basic_monomial<comp>> minimal_lcms;
        for (size_t idx = 0; idx < new_pairs.size(); ++idx)
        {
            bool eliminated = false;
            basic_monomial<comp> tmp;
            for (size_t midx = 0; midx < minimal_lcms.size(); ++midx)
            {
                if (is_divexact(tmp, new_lcms[idx], minimal_lcms[midx]))
                {
                    eliminated = true;
                    break;
                }
            }
            if (!eliminated)
            {
                minimal_new.push_back(new_pairs[idx]);
                minimal_lcms.push_back(new_lcms[idx]);
            }
        }
        new_pairs = std::move(minimal_new);

        // 4. 乘积准则: gcd(LM(G[k]), LM(h)) = 1 的对丢弃
        std::vector<__critical_pair> final_new;
        for (auto& p : new_pairs)
        {
            auto g = gcd(G[p.i].front().first, LM_h);
            if (g.deg() != 0)
                final_new.push_back(p);
        }
        new_pairs = std::move(final_new);

        // 5. 将 h 加入 G，新对加入 pairs
        G.push_back(h);
        sugar_vec.push_back(h_sugar);
        pairs.insert(pairs.end(), new_pairs.begin(), new_pairs.end());
    }

    // ===== M5: Buchberger 主循环 =====

    template<class Tc, class comp>
    inline void __interreduce(std::vector<polynomial_<Tc, comp>>& G)
    {
        // 步骤 1: 删除冗余
        size_t i = 0;
        while (i < G.size())
        {
            bool redundant = false;
            for (size_t j = 0; j < G.size(); ++j)
            {
                if (j == i) continue;
                basic_monomial<comp> tmp;
                if (is_divexact(tmp, G[i].front().first, G[j].front().first))
                {
                    redundant = true;
                    break;
                }
            }
            if (redundant)
                G.erase(G.begin() + i);
            else
                ++i;
        }

        // 步骤 2: 完全约化每个元素
        for (size_t i = 0; i < G.size(); ++i)
        {
            std::vector<polynomial_<Tc, comp>> others;
            others.reserve(G.size() - 1);
            for (size_t j = 0; j < G.size(); ++j)
                if (j != i) others.push_back(G[j]);
            G[i] = __normal_form(G[i], others);
            __make_monic(G[i]);
        }
    }

    template<class Tc, class comp>
    inline std::vector<polynomial_<Tc, comp>> __buchberger(
        std::vector<polynomial_<Tc, comp>> F)
    {
        // 1. 预处理: 过滤零多项式，首一化
        std::vector<polynomial_<Tc, comp>> G;
        for (auto& f : F)
        {
            if (!f.empty())
            {
                __make_monic(f);
                G.push_back(std::move(f));
            }
        }
        if (G.empty()) return G;
        if (G.size() == 1)
        {
            __interreduce(G);
            return G;
        }

        // 2. 初始化 sugar 度数
        std::vector<int64_t> sugar_vec;
        for (auto& g : G)
            sugar_vec.push_back(degree(g));

        // 3. 初始化临界对（逐个加入，复用 __update）
        std::vector<__critical_pair> pairs;
        std::vector<polynomial_<Tc, comp>> init_G;
        std::vector<int64_t> init_sugar;
        init_G.push_back(G[0]);
        init_sugar.push_back(sugar_vec[0]);
        for (size_t k = 1; k < G.size(); ++k)
            __update(init_G, pairs, G[k], init_sugar, sugar_vec[k]);
        G = std::move(init_G);
        sugar_vec = std::move(init_sugar);

        // 4. 主循环
        __cp_compare cmp;
        while (!pairs.empty())
        {
            auto min_it = std::min_element(pairs.begin(), pairs.end(), cmp);
            __critical_pair cp = *min_it;
            pairs.erase(min_it);

            auto h = __s_polynomial(G[cp.i], G[cp.j]);
            int64_t h_sugar = cp.sugar;
            h = __normal_form_with_sugar(h, G, sugar_vec, h_sugar);

            if (!h.empty())
            {
                __make_monic(h);
                __update(G, pairs, h, sugar_vec, h_sugar);
            }
        }

        // 5. 互约
        __interreduce(G);

        return G;
    }

    // ===== M6: 公开 API =====

    // QQ 入口
    template<class comp>
    inline std::vector<polynomial_<QQ, comp>> groebner_basis(
        const std::vector<polynomial_<QQ, comp>>& generators)
    {
        return __buchberger(generators);
    }

    // ZZ 入口（ZZ→QQ→ZZ 工作流）
    template<class comp>
    inline std::vector<polynomial_<ZZ, comp>> groebner_basis(
        const std::vector<polynomial_<ZZ, comp>>& F_zz)
    {
        // 1. ZZ → QQ
        std::vector<polynomial_<QQ, comp>> F_qq;
        F_qq.reserve(F_zz.size());
        for (auto& f : F_zz)
        {
            polynomial_<QQ, comp> fq(f.comp_ptr());
            poly_convert(f, fq);
            F_qq.push_back(std::move(fq));
        }

        // 2. 在 QQ 上计算
        auto G_qq = __buchberger(std::move(F_qq));

        // 3. QQ → ZZ（清分母 + 本原化）
        std::vector<polynomial_<ZZ, comp>> G_zz;
        G_zz.reserve(G_qq.size());
        for (auto& g : G_qq)
        {
            polynomial_<ZZ, comp> gz(g.comp_ptr());
            poly_convert(g, gz);

            // 本原化：计算所有系数的整数 GCD
            ZZ c = abs(gz.front().second);
            for (auto& term : gz.data())
            {
                c = gcd(c, abs(term.second));
                if (c == 1) break;
            }
            if (gz.front().second < 0)
                c = -c;  // 确保 lc > 0
            for (auto& term : gz.data())
                term.second /= c;

            G_zz.push_back(std::move(gz));
        }

        return G_zz;
    }

    // Ideal 类
    template<class Tc, class comp = grlex>
    class Ideal {
        std::vector<polynomial_<Tc, comp>> generators_;
        mutable std::vector<polynomial_<Tc, comp>> gb_cache_;
        mutable bool gb_computed_ = false;

    public:
        Ideal(std::vector<polynomial_<Tc, comp>> gens)
            : generators_(std::move(gens)) {}

        Ideal(std::initializer_list<polynomial_<Tc, comp>> gens)
            : generators_(gens) {}

        const std::vector<polynomial_<Tc, comp>>& generators() const {
            return generators_;
        }

        const std::vector<polynomial_<Tc, comp>>& groebner_basis() const {
            if (!gb_computed_) {
                gb_cache_ = clpoly::groebner_basis(generators_);
                gb_computed_ = true;
            }
            return gb_cache_;
        }

        bool contains(const polynomial_<Tc, comp>& f) const {
            auto nf = normal_form(f, groebner_basis());
            return nf.empty();
        }

        polynomial_<Tc, comp> reduce(const polynomial_<Tc, comp>& f) const {
            return normal_form(f, groebner_basis());
        }
    };

} // namespace clpoly
#endif
