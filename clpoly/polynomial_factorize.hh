/**
 * @file polynomial_factorize.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 多项式因式分解（公开 API）
 *
 * 本文件是因式分解模块的唯一公开入口。内部实现分布在：
 *   - polynomial_factorize_zp.hh    : Zp[x] 辅助 + M1 (squarefree/DDF/EDF)
 *   - polynomial_factorize_univar.hh: ZZ[x] 辅助 + M2 Hensel + M3 Zassenhaus
 *   - polynomial_factorize_wang.hh  : 多变量 Wang 算法
 */
#ifndef CLPOLY_POLYNOMIAL_FACTORIZE_HH
#define CLPOLY_POLYNOMIAL_FACTORIZE_HH

#include <clpoly/polynomial_factorize_wang.hh>

namespace clpoly{

    // ================================================================
    // §8.5 factorize: ZZ[x] lex 特化
    // ================================================================

    template<class var_order>
    factorization<polynomial_<ZZ,lex_<var_order>>>
    factorize(const polynomial_<ZZ,lex_<var_order>>& F)
    {
        using Poly = polynomial_<ZZ,lex_<var_order>>;
        factorization<Poly> result;
        result.content = ZZ(1);

        // 零多项式
        if (F.empty())
        {
            result.content = ZZ(0);
            return result;
        }

        // 常数多项式
        if (is_number(F))
        {
            result.content = F.front().second;
            return result;
        }

        // 多变量 → dispatch 到 __factor_multivar
        auto vars = get_variables(F);
        if (vars.size() > 1)
            return __factor_multivar(F);

        variable var = vars.front().first;

        // 转 upolynomial
        upolynomial_<ZZ> uf;
        poly_convert(F, uf);

        // 提取内容，本原化
        auto [ct, uf_prim] = __upoly_primitive(std::move(uf));
        result.content = ct;

        // 线性多项式
        if (get_deg(uf_prim) <= 1)
        {
            auto p = __upoly_to_poly(uf_prim, var, F.comp_ptr());
            result.factors.push_back({std::move(p), 1});
            return result;
        }

        // 转回 polynomial 做无平方分解
        auto poly_prim = __upoly_to_poly(uf_prim, var, F.comp_ptr());
        auto sqf = squarefreefactorize(poly_prim);

        for (auto& [sqf_factor, mult] : sqf)
        {
            if (is_number(sqf_factor))
            {
                result.content *= sqf_factor.front().second;
                continue;
            }

            // 转 upolynomial
            upolynomial_<ZZ> usqf;
            poly_convert(sqf_factor, usqf);

            if (get_deg(usqf) <= 1)
            {
                // 线性因子直接加入
                auto p = __upoly_to_poly(usqf, var, F.comp_ptr());
                result.factors.push_back({std::move(p), mult});
                continue;
            }

            // 对 deg>=2 的无平方因子做不可约分解
            auto irr_factors = __factor_squarefree_primitive_ZZ(usqf);
            for (auto& irr : irr_factors)
            {
                // 确保 lc > 0
                if (irr.front().second < ZZ(0))
                {
                    for (auto& term : irr)
                        term.second = -term.second;
                }
                auto p = __upoly_to_poly(irr, var, F.comp_ptr());
                result.factors.push_back({std::move(p), mult});
            }
        }

        // 按 degree 排序
        std::sort(result.factors.begin(), result.factors.end(),
            [](const auto& a, const auto& b) {
                return degree(a.first) < degree(b.first);
            });

#ifndef NDEBUG
        {
            Poly check(F.comp_ptr());
            check.data().push_back({{}, result.content});
            check.normalization();
            for (const auto& [fi, ei] : result.factors)
                for (uint64_t e = 0; e < ei; ++e)
                {
                    check = check * fi;
                    check.normalization();
                }
            assert(check == F);
        }
#endif

        return result;
    }

    // ================================================================
    // §8.6 factorize: ZZ[x] 通用 comp 包装
    // ================================================================

    template<class comp>
    factorization<polynomial_<ZZ,comp>>
    factorize(const polynomial_<ZZ,comp>& F)
    {
        using Poly = polynomial_<ZZ,comp>;
        // 转 lex
        polynomial_<ZZ,lex> F_lex;
        poly_convert(F, F_lex);
        auto result_lex = factorize(F_lex);

        // 转回原 comp
        factorization<Poly> result;
        result.content = std::move(result_lex.content);
        for (auto& [fac, mult] : result_lex.factors)
        {
            Poly fac_comp(F.comp_ptr());
            poly_convert(fac, fac_comp);
            result.factors.push_back({std::move(fac_comp), mult});
        }
        return result;
    }

    // ================================================================
    // §8.7 factorize: QQ[x]
    // ================================================================

    template<class comp>
    factorization<polynomial_<QQ,comp>>
    factorize(const polynomial_<QQ,comp>& F)
    {
        using PolyQQ = polynomial_<QQ,comp>;
        factorization<PolyQQ> result;
        result.content = QQ(1);

        if (F.empty())
        {
            result.content = QQ(0);
            return result;
        }

        if (is_number(F))
        {
            result.content = F.front().second;
            return result;
        }

        // poly_convert(QQ→ZZ) 将系数乘以 LCD，得 F_zz = lcd · F
        polynomial_<ZZ,comp> F_zz(F.comp_ptr());
        poly_convert(F, F_zz);

        // ZZ 分解: F_zz = content_zz · ∏ fi^ei
        auto result_zz = factorize(F_zz);

        // 计算 LCD
        ZZ lcd(1);
        for (auto& term : F)
            lcd = lcm(lcd, term.second.get_den());

        // F = F_zz / lcd = (content_zz / lcd) · ∏ fi^ei
        // 令 fi = lc(fi) · fi_monic，则:
        //   F = (content_zz / lcd) · ∏ lc(fi)^ei · ∏ fi_monic^ei
        // QQ content = (content_zz / lcd) · ∏ lc(fi)^ei
        result.content = QQ(result_zz.content, lcd);

        for (auto& [fac_zz, mult] : result_zz.factors)
        {
            // fi / lc(fi) → 首一 QQ 因子
            PolyQQ fac_qq(F.comp_ptr());
            ZZ lc_val = fac_zz.front().second;
            for (auto& term : fac_zz)
            {
                basic_monomial<comp> m(F.comp_ptr());
                m = term.first;
                fac_qq.push_back({std::move(m), QQ(term.second, lc_val)});
            }

            // lc(fi)^ei 吸收进 content
            result.content *= QQ(pow(lc_val, (int64_t)mult), ZZ(1));

            result.factors.push_back({std::move(fac_qq), mult});
        }

#ifndef NDEBUG
        {
            PolyQQ check(F.comp_ptr());
            check.data().push_back({{}, result.content});
            check.normalization();
            for (const auto& [fi, ei] : result.factors)
                for (uint64_t e = 0; e < ei; ++e)
                {
                    check = check * fi;
                    check.normalization();
                }
            assert(check == F);
        }
#endif

        return result;
    }
    // ================================================================
    // factorbasis: 用不可约分解构造互素基
    // 接口与 squarefreebasis 一致
    // 返回: (基多项式列表, 每个基元素的来源 {(输入下标, 重数), ...})
    // ================================================================

    // 字典序核心实现
    template<class var_order>
    std::pair<std::vector<polynomial_<ZZ,lex_<var_order>>>,
              std::vector<std::vector<std::pair<uint64_t,uint64_t>>>>
    factorbasis(const std::vector<polynomial_<ZZ,lex_<var_order>>>& F)
    {
        using Poly = polynomial_<ZZ,lex_<var_order>>;
        std::vector<Poly> lst;
        std::vector<std::vector<std::pair<uint64_t,uint64_t>>> I;

        for (size_t i = 0; i < F.size(); ++i)
        {
            auto fac = factorize(F[i]);

            for (auto& [factor, mult] : fac.factors)
            {
                // factorize 保证 ZZ 因子非常数、lc > 0 且本原，直接用 == 比较
                assert(!is_number(factor));
                bool found = false;
                for (size_t j = 0; j < lst.size(); ++j)
                {
                    if (lst[j] == factor)
                    {
                        I[j].push_back({i, mult});
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    lst.push_back(std::move(factor));
                    I.push_back({{i, mult}});
                }
            }
        }

        return {std::move(lst), std::move(I)};
    }

    // 一般序包装：转 lex → 调核心 → 转回
    template<class comp>
    std::pair<std::vector<polynomial_<ZZ,comp>>,
              std::vector<std::vector<std::pair<uint64_t,uint64_t>>>>
    factorbasis(const std::vector<polynomial_<ZZ,comp>>& F)
    {
        if (F.empty())
            return {{}, {}};

        // 转字典序
        std::vector<polynomial_<ZZ,lex>> lst;
        polynomial_<ZZ,lex> p;
        lst.reserve(F.size());
        for (const auto& i : F)
        {
            poly_convert(i, p);
            lst.push_back(std::move(p));
        }

        auto [basis_lex, info] = factorbasis(lst);

        // 转回原序
        std::vector<polynomial_<ZZ,comp>> basis;
        basis.reserve(basis_lex.size());
        polynomial_<ZZ,comp> p1(F.front().comp_ptr());
        for (auto& i : basis_lex)
        {
            poly_convert(i, p1);
            basis.push_back(std::move(p1));
        }

        return {std::move(basis), std::move(info)};
    }

} // namespace clpoly

#endif // CLPOLY_POLYNOMIAL_FACTORIZE_HH
