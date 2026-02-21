/**
 * @file polynomial_factorize_wang.hh
 * @brief 多变量因式分解：Wang 算法辅助、多变量 Hensel 提升、
 *        Wang 核心 (__wang_core)、多变量入口 (__factor_multivar)。
 *
 * 本文件可独立包含。供 polynomial_factorize.hh 使用。
 *
 * 注意：__factor_multivar 内部递归调用 factorize(polynomial_<ZZ,lex_>)，
 *       该函数定义在 polynomial_factorize.hh 中，通过 C++ 两阶段查找
 *       在实例化时解析（依赖调用，合法）。
 */
#ifndef CLPOLY_POLYNOMIAL_FACTORIZE_WANG_HH
#define CLPOLY_POLYNOMIAL_FACTORIZE_WANG_HH

#include <clpoly/polynomial_factorize_univar.hh>
#include <map>
#include <set>

namespace clpoly{

    // Taylor 系数提取: f 在 xk = alpha_k 处的第 j 阶 Taylor 系数
    // 将 f 写成 Σ cⱼ·(xk - alpha_k)^j，返回 cⱼ
    template<class var_order>
    polynomial_<ZZ, lex_<var_order>> __taylor_coeff(
        const polynomial_<ZZ, lex_<var_order>>& f,
        const variable& xk, const ZZ& alpha_k, int j)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        Poly g = f;
        // 构造 (xk - alpha_k) 作为 polynomial
        Poly divisor(f.comp_ptr());
        {
            basic_monomial<lex_<var_order>> m_xk(f.comp_ptr());
            m_xk.push_back({xk, 1});
            divisor.push_back({m_xk, ZZ(1)});
            if (alpha_k != 0)
            {
                basic_monomial<lex_<var_order>> m_const(f.comp_ptr());
                divisor.push_back({m_const, -alpha_k});
            }
            divisor.normalization();
        }
        // 反复精确除以 (xk - alpha_k) j 次
        for (int t = 0; t < j; ++t)
        {
            Poly q(f.comp_ptr()), r(f.comp_ptr());
            pair_vec_div(q.data(), r.data(), g.data(), divisor.data(), f.comp());
            g = std::move(q);
        }
        // 代入 xk = alpha_k
        return assign(g, xk, alpha_k);
    }
    // ================================================================
    // §M5 多变量因式分解辅助
    // ================================================================

    // §4.2 选取值点: 返回 α = {x₂→α₂,...,xₙ→αₙ}
    // 满足 (a) f(x₁,α) 无平方 (b) lc(f,x₁)(α)≠0 (d) lc 因子求值两两互素
    template<class var_order>
    std::map<variable, ZZ>
    __select_eval_point(
        const polynomial_<ZZ, lex_<var_order>>& f,
        const variable& main_var,
        int skip = 0)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;

        // 收集除主变量外的变量
        auto all_vars = get_variables(f);
        std::vector<variable> vars;
        for (auto& [v, d] : all_vars)
            if (v != main_var)
                vars.push_back(v);
        assert(!vars.empty());
        size_t n = vars.size();

        // L = lc(f, x₁)
        Poly L = leadcoeff(f, main_var);

        // 预处理: 若 L 非常数，分解 L 用于条件 (d)
        std::vector<Poly> lc_irr_factors;
        if (!is_number(L))
        {
            auto lc_fac = factorize(L);
            for (auto& [lj, ej] : lc_fac.factors)
                lc_irr_factors.push_back(lj);
        }

        // 枚举候选点
        for (int bound = 0; ; ++bound)
        {
            // 生成所有 n 维向量, 各分量在 [-bound, bound], max(|αᵢ|) == bound
            // 总数 (2*bound+1)^n - (2*(bound-1)+1)^n (排除内层)
            int range = 2 * bound + 1;
            int total = 1;
            for (size_t i = 0; i < n; ++i) total *= range;

            for (int idx = 0; idx < total; ++idx)
            {
                // 解码 idx 为各分量值
                std::map<variable, ZZ> alpha;
                int tmp = idx;
                int max_abs = 0;
                for (size_t i = 0; i < n; ++i)
                {
                    int val = (tmp % range) - bound;
                    tmp /= range;
                    alpha[vars[i]] = ZZ(val);
                    if (std::abs(val) > max_abs) max_abs = std::abs(val);
                }
                // 跳过已在更小 bound 检查过的点
                if (max_abs < bound) continue;

                // 条件 (b): lc 非零
                ZZ delta;
                if (is_number(L))
                    delta = L.front().second;
                else
                {
                    auto L_eval = assign(L, alpha);
                    if (L_eval.empty()) continue;
                    if (!is_number(L_eval)) continue;
                    delta = L_eval.front().second;
                }
                if (delta == 0) continue;

                // 条件 (a): f(x₁, α) 无平方
                // 注: f0 是单变量多项式, 用 upolynomial GCD 避免多变量 GCD 的开销
                auto f0 = assign(f, alpha);
                {
                    upolynomial_<ZZ> f0_u;
                    poly_convert(f0, f0_u);
                    auto f0_d = derivative(f0_u);
                    if (!f0_d.empty())
                    {
                        auto g = polynomial_GCD(f0_u, f0_d);
                        if (!g.empty() && get_deg(g) > 0) continue;
                    }
                }

                // 条件 (d): lc 因子求值两两互素
                if (!lc_irr_factors.empty())
                {
                    std::vector<ZZ> vals;
                    bool bad_point = false;
                    for (auto& lj : lc_irr_factors)
                    {
                        auto lj_eval = assign(lj, alpha);
                        if (lj_eval.empty()) { bad_point = true; break; }
                        ZZ v = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
                        if (v == 0) { bad_point = true; break; }
                        vals.push_back(abs(v));
                    }
                    if (bad_point) continue;
                    // 条件 (d'): |lⱼ(α)| ≥ 2, 避免退化赋值
                    for (auto& v : vals)
                        if (v <= ZZ(1)) { bad_point = true; break; }
                    if (bad_point) continue;
                    // 两两互素
                    bool pairwise_coprime = true;
                    for (size_t i = 0; i < vals.size() && pairwise_coprime; ++i)
                        for (size_t j = i + 1; j < vals.size() && pairwise_coprime; ++j)
                            if (gcd(vals[i], vals[j]) != ZZ(1))
                                pairwise_coprime = false;
                    if (!pairwise_coprime) continue;
                }

                if (skip > 0) { --skip; continue; }
                return alpha;
            }
        }
        // 不可达
        return {};
    }

    // §5.4 首项系数校正结果
    template<class var_order>
    struct __wang_lc_result {
        bool success;
        polynomial_<ZZ, lex_<var_order>> f_scaled;
        std::vector<polynomial_<ZZ, lex_<var_order>>> lc_assignments;  // σᵢ
        std::vector<polynomial_<ZZ, lex_<var_order>>> lc_targets;      // τᵢ
        std::vector<upolynomial_<ZZ>> scaled_factors;                  // vᵢ
    };

    // §5.2 首项系数校正
    template<class var_order>
    __wang_lc_result<var_order> __wang_leading_coeff(
        const polynomial_<ZZ, lex_<var_order>>& f,
        const std::vector<upolynomial_<ZZ>>& univar_factors,
        const std::map<variable, ZZ>& eval_point,
        const variable& main_var,
        const ZZ& uni_content = ZZ(1))
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        auto comp_ptr = f.comp_ptr();

        auto make_const = [&](const ZZ& val) -> Poly {
            Poly p(comp_ptr);
            if (val != 0)
                p.push_back({basic_monomial<lex_<var_order>>(comp_ptr), val});
            return p;
        };

        size_t r = univar_factors.size();
        __wang_lc_result<var_order> result;
        result.success = false;

        // Step 1: L = lc(f, x₁), δ = L(α)
        Poly L = leadcoeff(f, main_var);
        ZZ delta;
        if (is_number(L))
            delta = L.front().second;
        else
        {
            auto L_eval = assign(L, eval_point);
            if (L_eval.empty() || !is_number(L_eval)) return result;
            delta = L_eval.front().second;
        }
        if (delta == 0) return result;

        // σᵢ 初始化
        std::vector<Poly> sigma(r, make_const(ZZ(1)));

        if (!is_number(L))
        {
            // Step 2: 递归分解 L
            auto lc_fac = factorize(L);
            ZZ gamma = lc_fac.content;

            // Step 3: 分配 lc 因子
            // 收集 (lⱼ, eⱼ) 并按 eⱼ 降序排序
            auto& lc_factors = lc_fac.factors;
            std::sort(lc_factors.begin(), lc_factors.end(),
                [](const auto& a, const auto& b) { return a.second > b.second; });

            // wᵢ = lc(uᵢ) * cs — 数值跟踪 (cs = content(f₀))
            std::vector<ZZ> w(r);
            for (size_t i = 0; i < r; ++i)
                w[i] = univar_factors[i].front().second * uni_content;

            // P2: 验证恒等式 δ = γ·∏lⱼ(α)^eⱼ = cs·∏lc(uᵢ)
            {
                ZZ prod_z = gamma;
                for (auto& [lj2, ej2] : lc_factors)
                {
                    auto lj2_eval = assign(lj2, eval_point);
                    ZZ v = is_number(lj2_eval) ? lj2_eval.front().second : ZZ(0);
                    for (uint64_t e = 0; e < ej2; ++e)
                        prod_z *= v;
                }
                ZZ prod_lc(uni_content);
                for (size_t i = 0; i < r; ++i)
                    prod_lc *= univar_factors[i].front().second;
                if (abs(prod_z) != abs(prod_lc))
                    return result;  // 数据不一致，换点
            }

            // GCD 匹配 LC 因子 (SymPy dmp_zz_wang_lead_coeffs 风格)
            // 对每个 LC 因子 lⱼ^eⱼ, 找 gcd(|w[i]|, |lⱼ(α)|^eⱼ) 最大的唯一 i
            for (auto& [lj, ej] : lc_factors)
            {
                auto lj_eval = assign(lj, eval_point);
                ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
                if (lj_val == 0 || abs(lj_val) == ZZ(1)) return result;

                ZZ lj_pow(1);
                for (uint64_t e = 0; e < ej; ++e)
                    lj_pow = lj_pow * abs(lj_val);

                // 找唯一最佳匹配
                size_t best_i = r;
                ZZ best_g(0);
                bool ambiguous = false;
                for (size_t i = 0; i < r; ++i)
                {
                    ZZ g = gcd(abs(w[i]), lj_pow);
                    if (g > best_g)
                    {
                        best_g = g;
                        best_i = i;
                        ambiguous = false;
                    }
                    else if (g == best_g && g > ZZ(1))
                    {
                        ambiguous = true;
                    }
                }
                if (best_i >= r || ambiguous || best_g <= ZZ(1))
                    return result;

                // σ[best_i] *= lⱼ^eⱼ
                Poly lj_power = lj;
                for (uint64_t e = 1; e < ej; ++e)
                {
                    lj_power = lj_power * lj;
                    lj_power.normalization();
                }
                sigma[best_i] = sigma[best_i] * lj_power;
                sigma[best_i].normalization();

                // 从 w 中移除已匹配的部分
                w[best_i] /= best_g;
            }

            // 吸收整数内容 γ 到 σ₀
            sigma[0] = sigma[0] * make_const(gamma);
            sigma[0].normalization();
        }
        else
        {
            // L 是常数 → σ₁ = L, 其余 = 1
            sigma[0] = L;
        }

        // Step 4: 缩放
        // f_scaled = δ^(r-1) * f
        ZZ delta_pow(1);
        for (size_t i = 1; i < r; ++i)
            delta_pow *= delta;
        result.f_scaled = f * make_const(delta_pow);
        result.f_scaled.normalization();

        // vᵢ = δ * (uᵢ / lc(uᵢ))
        result.scaled_factors.resize(r);
        for (size_t i = 0; i < r; ++i)
        {
            ZZ lc_ui = univar_factors[i].front().second;
            // 首一化: uᵢ / lc(uᵢ) 在 Z[x] 上不精确 → 用 δ/lc(uᵢ) * uᵢ
            // vᵢ = δ * ūᵢ = (δ / lc(uᵢ)) * uᵢ
            // 注: δ / lc(uᵢ) 不一定整除，所以分两步:
            // vᵢ 的每项系数 = δ * coeff / lc(uᵢ)
            upolynomial_<ZZ> vi;
            vi.reserve(univar_factors[i].size());
            for (auto& term : univar_factors[i])
            {
                assert(delta * term.second % lc_ui == ZZ(0));
                ZZ new_coeff = delta * term.second / lc_ui;
                if (new_coeff != 0)
                    vi.push_back({term.first, new_coeff});
            }
            vi.normalization();
            result.scaled_factors[i] = std::move(vi);
        }

        // Step 5: τᵢ = (δ / σᵢ(α)) * σᵢ
        result.lc_assignments = sigma;
        result.lc_targets.resize(r);
        for (size_t i = 0; i < r; ++i)
        {
            auto sigma_eval = assign(sigma[i], eval_point);
            ZZ sigma_val = is_number(sigma_eval) ? sigma_eval.front().second : ZZ(0);
            if (sigma_val == 0) return result;  // 不应发生
            assert(delta % sigma_val == ZZ(0));
            ZZ scale = delta / sigma_val;
            result.lc_targets[i] = sigma[i] * make_const(scale);
            result.lc_targets[i].normalization();
        }

        result.success = true;
        return result;
    }

    // ================================================================
    // §6 多变量 Hensel 提升
    // ================================================================

    // §6.0 辅助: 多变量多项式对 main_var 做伪除法
    // 将 f ∈ Z[x₁,...,xₙ] 视为 R[x₁] 其中 R = Z[x₂,...,xₙ]
    // 对 g ∈ Z[x₁] (单变量) 做伪除法:  lc(g)^ki · f = q·g + rem
    // 返回 (rem, lc(g)^ki), 其中 ki = max(deg_x₁(f) - deg_x₁(g) + 1, 0)
    template<class var_order>
    std::pair<polynomial_<ZZ, lex_<var_order>>, ZZ>
    __pseudo_remainder_x1(
        const polynomial_<ZZ, lex_<var_order>>& f,
        const upolynomial_<ZZ>& g,
        const variable& main_var)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        using BMono = basic_monomial<lex_<var_order>>;
        auto comp_ptr = f.comp_ptr();

        int64_t m = get_deg(g);
        ZZ lc_g = g.front().second;
        ZZ scale(1);

        // 按 main_var 的次数分组: coeffs[d] ∈ Z[其他变量]
        std::map<int64_t, Poly> coeffs;
        for (auto& term : f)
        {
            int64_t x1_deg = 0;
            BMono rest_mono(comp_ptr);
            for (auto& [var, deg] : term.first)
            {
                if (var == main_var)
                    x1_deg = deg;
                else
                    rest_mono.push_back({var, (uint64_t)deg});
            }
            rest_mono.normalization();
            if (coeffs.find(x1_deg) == coeffs.end())
                coeffs[x1_deg] = Poly(comp_ptr);
            coeffs[x1_deg].push_back({rest_mono, term.second});
        }
        for (auto& [d, p] : coeffs)
            p.normalization();

        // 伪除法主循环
        while (true)
        {
            // 找最高非空 x₁ 次数
            int64_t d = -1;
            for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it)
                if (!it->second.empty()) { d = it->first; break; }
            if (d < m) break;

            Poly c_d = coeffs[d];

            // 缩放所有系数乘以 lc(g)
            for (auto& [deg, p] : coeffs)
            {
                for (auto& term : p.data())
                    term.second *= lc_g;
                p.normalization();
            }
            scale *= lc_g;

            // 减去 c_d · x₁^(d-m) · g
            for (auto& g_term : g)
            {
                int64_t g_deg = g_term.first.deg();
                ZZ g_coeff = g_term.second;
                int64_t target_deg = d - m + g_deg;
                // coeffs[target_deg] -= c_d * g_coeff
                Poly subtrahend(comp_ptr);
                for (auto& [mono, coeff] : c_d)
                    subtrahend.push_back({mono, coeff * g_coeff});
                subtrahend.normalization();

                if (coeffs.find(target_deg) == coeffs.end())
                    coeffs[target_deg] = Poly(comp_ptr);
                coeffs[target_deg] = coeffs[target_deg] - subtrahend;
                coeffs[target_deg].normalization();
            }
        }

        // 重组余式
        Poly rem(comp_ptr);
        for (auto& [d, p] : coeffs)
        {
            for (auto& [mono, coeff] : p)
            {
                if (coeff == ZZ(0)) continue;
                BMono full_mono(comp_ptr);
                if (d > 0)
                    full_mono.push_back({main_var, (uint64_t)d});
                for (auto& [var, deg] : mono)
                    full_mono.push_back({var, (uint64_t)deg});
                full_mono.normalization();
                rem.push_back({full_mono, coeff});
            }
        }
        rem.normalization();

        return {rem, scale};
    }


    // §6.2 递归多变量 Diophantine 求解器 (GCL Algorithm 6.3)
    // 求解: h = Σ δᵢ · Ĝᵢ, 其中 Ĝᵢ = ∏_{j≠i} G[j]
    // eval_vars: 已提升的变量及其求值点 [(x₂,α₂),...,(x_{k-1},α_{k-1})]
    // depth: 当前递归深度 (从 0 开始剥离 eval_vars)
    template<class var_order>
    std::vector<polynomial_<ZZ, lex_<var_order>>>
    __multivar_diophantine(
        const polynomial_<ZZ, lex_<var_order>>& h,
        const std::vector<polynomial_<ZZ, lex_<var_order>>>& G,
        const std::vector<upolynomial_<ZZ>>& bezout_s,
        const ZZ& bezout_denom,
        const std::vector<upolynomial_<ZZ>>& v_factors,
        const variable& main_var,
        const std::vector<std::pair<variable, ZZ>>& eval_vars,
        const ZZ& pa,
        size_t depth)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        auto comp_ptr = h.comp_ptr();
        size_t r = G.size();

        if (depth >= eval_vars.size())
        {
            // 基本情形: h 是单变量 (只含 main_var), 用 Bézout 链
            std::vector<Poly> result(r, Poly(comp_ptr));
            upolynomial_<ZZ> h_upoly;
            poly_convert(h, h_upoly);

            for (size_t i = 0; i < r; ++i)
            {
                auto si_h = bezout_s[i] * h_upoly;
                si_h.normalization();
                if (si_h.empty()) continue;

                // 伪除法: lc(vi)^k * (si*h) = q * vi + rem, k = deg(si_h)-deg(vi)+1
                int64_t k = std::max(get_deg(si_h) - get_deg(v_factors[i]) + 1, (int64_t)0);

                upolynomial_<ZZ> rem_upoly;
                upoly_prem(rem_upoly, si_h, v_factors[i], main_var);

                ZZ lc_vi = v_factors[i].front().second;
                ZZ lc_pow(1);
                for (int64_t e = 0; e < k; ++e)
                    lc_pow *= lc_vi;

                // 总缩放因子 = bezout_denom * lc_pow
                ZZ divisor = bezout_denom * lc_pow;
                upolynomial_<ZZ> delta_i_upoly;

                if (pa > ZZ(0))
                {
                    // 模逆路径: divisor⁻¹ mod pa, 然后对称模约化
                    ZZ divisor_inv;
                    bool ok = ZZ::invert(divisor_inv, divisor, pa);
                    assert(ok);  // gcd(p, divisor) = 1 由选素数保证
                    for (auto& term : rem_upoly)
                    {
                        ZZ coeff = __symmetric_mod(term.second * divisor_inv, pa);
                        if (coeff != ZZ(0))
                            delta_i_upoly.push_back({term.first, coeff});
                    }
                }
                else
                {
                    // 精确除法路径 (denom == ±1)
                    for (auto& term : rem_upoly)
                    {
                        ZZ coeff = term.second / divisor;
                        if (coeff != ZZ(0))
                            delta_i_upoly.push_back({term.first, coeff});
                    }
                }

                poly_convert(delta_i_upoly, result[i], main_var);
            }
            return result;
        }

        // 递归情形: 剥离 eval_vars[depth] 对应的变量
        auto [xk, alpha_k] = eval_vars[depth];

        // Step 1: 求值 h → h_bar (少一个变量)
        std::map<variable, ZZ> sub = {{xk, alpha_k}};
        Poly h_bar = assign(h, sub);

        // Step 2: 递归求解 (少一个变量)
        // 传入求值后的 G (但注意：递归传同样的 G，因为内层只用到 Bézout 链)
        auto delta = __multivar_diophantine(
            h_bar, G, bezout_s, bezout_denom, v_factors,
            main_var, eval_vars, pa, depth + 1);

        // Step 3: 计算误差 e = h - Σ δᵢ · Ĝᵢ
        // Ĝᵢ = ∏_{j≠i} G[j] (多变量乘积)
        Poly correction(comp_ptr);
        for (size_t i = 0; i < r; ++i)
        {
            if (delta[i].empty()) continue;
            // 计算 Ĝᵢ = ∏_{j≠i} G[j]
            Poly G_hat_i(comp_ptr);
            bool first = true;
            for (size_t j = 0; j < r; ++j)
            {
                if (j == i) continue;
                if (first)
                {
                    G_hat_i = G[j];
                    first = false;
                }
                else
                {
                    G_hat_i = G_hat_i * G[j];
                    G_hat_i.normalization();
                }
            }
            auto term = delta[i] * G_hat_i;
            term.normalization();
            correction = correction + term;
            correction.normalization();
        }
        Poly error = h - correction;
        error.normalization();

        if (error.empty()) return delta;

        // Step 4: 对 error 按 (xk - αk) 的 Taylor 展开, 逐阶求解
        int dk = degree(error, xk);
        Poly xk_minus_alpha(comp_ptr);
        {
            basic_monomial<lex_<var_order>> m_xk(comp_ptr);
            m_xk.push_back({xk, 1});
            xk_minus_alpha.push_back({m_xk, ZZ(1)});
            if (alpha_k != 0)
            {
                basic_monomial<lex_<var_order>> m0(comp_ptr);
                xk_minus_alpha.push_back({m0, -alpha_k});
            }
            xk_minus_alpha.normalization();
        }
        Poly xk_pow = xk_minus_alpha;  // (xk - αk)^1

        for (int j = 1; j <= dk; ++j)
        {
            auto cj = __taylor_coeff(error, xk, alpha_k, j);
            if (cj.empty())
            {
                if (j < dk)
                {
                    xk_pow = xk_pow * xk_minus_alpha;
                    xk_pow.normalization();
                }
                continue;
            }

            // 递归求解 cj (cj 已不含 xk)
            auto tau = __multivar_diophantine(
                cj, G, bezout_s, bezout_denom, v_factors,
                main_var, eval_vars, pa, depth + 1);

            for (size_t i = 0; i < r; ++i)
            {
                if (tau[i].empty()) continue;
                auto update = tau[i] * xk_pow;
                update.normalization();
                delta[i] = delta[i] + update;
                delta[i].normalization();
            }

            // GCL §6.3: 更新误差 (多变量余因子产生交叉项, 影响后续阶系数)
            // error -= Σ τᵢ · (xk-αk)^j · Ĝᵢ
            for (size_t i = 0; i < r; ++i)
            {
                if (tau[i].empty()) continue;
                Poly G_hat_i(comp_ptr);
                bool first = true;
                for (size_t l = 0; l < r; ++l)
                {
                    if (l == i) continue;
                    if (first) { G_hat_i = G[l]; first = false; }
                    else { G_hat_i = G_hat_i * G[l]; G_hat_i.normalization(); }
                }
                auto sub = tau[i] * xk_pow * G_hat_i;
                sub.normalization();
                error = error - sub;
                error.normalization();
            }
            if (error.empty()) break;

            if (j < dk)
            {
                xk_pow = xk_pow * xk_minus_alpha;
                xk_pow.normalization();
            }
        }

        return delta;
    }

    // §6.3 LC 校正: 将 G[i] 关于 main_var 的首项系数替换为 lc_target
    template<class var_order>
    void __hensel_lc_correct(
        polynomial_<ZZ, lex_<var_order>>& Gi,
        const polynomial_<ZZ, lex_<var_order>>& lc_target,
        const variable& main_var)
    {
        if (Gi.empty()) return;
        auto lc_current = leadcoeff(Gi, main_var);
        auto diff = lc_target - lc_current;
        diff.normalization();
        if (diff.empty()) return;

        int64_t d = degree(Gi, main_var);
        // 构造 main_var^d
        polynomial_<ZZ, lex_<var_order>> x1d(Gi.comp_ptr());
        basic_monomial<lex_<var_order>> m_x1d(Gi.comp_ptr());
        m_x1d.push_back({main_var, (uint64_t)d});
        x1d.push_back({m_x1d, ZZ(1)});

        auto correction = diff * x1d;
        correction.normalization();
        Gi = Gi + correction;
        Gi.normalization();
    }

    // §6.3 单变量 Hensel 提升步
    // Pre:  ∏ Gᵢ|_{xk=αk} = f_curr|_{xk=αk} (LC 校正后)
    //       Σ sᵢ·V̂ᵢ = bezout_denom, V̂ᵢ = ∏_{j≠i} vⱼ
    // Post: Gᵢ 扩展到包含 xk, ∏ Gᵢ = f_curr
    // prev_eval: 前序已提升变量的求值点 {x₂→α₂,...,x_{k-1}→α_{k-1}}
    template<class var_order>
    void __hensel_lift_one_var(
        const polynomial_<ZZ, lex_<var_order>>& f_curr,
        std::vector<polynomial_<ZZ, lex_<var_order>>>& G,
        const std::vector<upolynomial_<ZZ>>& bezout_s,
        const ZZ& bezout_denom,
        const ZZ& pa,
        const std::vector<upolynomial_<ZZ>>& v_factors,
        const ZZ& delta,
        const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_tau,
        const variable& main_var,
        const variable& xk, const ZZ& alpha_k, int dk,
        const std::map<variable, ZZ>& prev_eval)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        auto comp_ptr = f_curr.comp_ptr();
        size_t r = G.size();

        // 构造 (xk - αk)
        Poly xk_minus_alpha(comp_ptr);
        {
            basic_monomial<lex_<var_order>> m_xk(comp_ptr);
            m_xk.push_back({xk, 1});
            xk_minus_alpha.push_back({m_xk, ZZ(1)});
            if (alpha_k != 0)
            {
                basic_monomial<lex_<var_order>> m0(comp_ptr);
                xk_minus_alpha.push_back({m0, -alpha_k});
            }
            xk_minus_alpha.normalization();
        }

        // Step E (pre-step): LC 校正
        for (size_t i = 0; i < r; ++i)
            __hensel_lc_correct(G[i], lc_tau[i], main_var);

        // G_base = G|_{xk=αk}: Diophantine 求解器用的因子
        // 方程 ej = Σ δᵢ · ∏_{l≠i} Gₗ|_{xk=αk} 中需要去掉 xk 方向的贡献
        // (xk 修正 ∝ (xk-αk)^m 在 xk=αk 处为零, 故 G_base 在循环中不变)
        std::vector<Poly> G_base(r, Poly(comp_ptr));
        {
            std::map<variable, ZZ> xk_sub = {{xk, alpha_k}};
            for (size_t i = 0; i < r; ++i)
                G_base[i] = assign(G[i], xk_sub);
        }

        // 构建 eval_vars: 已提升变量的求值点列表 (循环中不变)
        std::vector<std::pair<variable, ZZ>> eval_vars;
        for (auto& [v, val] : prev_eval)
            eval_vars.push_back({v, val});

        // (xk - αk)^j, 从 j=1 开始
        Poly xk_pow = xk_minus_alpha;  // (xk - αk)^1

        // 缓存 ∏Gᵢ, 仅在 G 实际更新后重算
        Poly prod = G[0];
        for (size_t i = 1; i < r; ++i)
        {
            prod = prod * G[i];
            prod.normalization();
        }
        bool prod_dirty = false;

        for (int j = 1; j <= dk; ++j)
        {
            // Step A: e = f_curr - ∏Gᵢ
            if (prod_dirty)
            {
                prod = G[0];
                for (size_t i = 1; i < r; ++i)
                {
                    prod = prod * G[i];
                    prod.normalization();
                }
                prod_dirty = false;
            }
            Poly e = f_curr - prod;
            e.normalization();

            if (e.empty()) break;  // 提前终止

            // Step B: eⱼ = __taylor_coeff(e, xk, αk, j)
            auto ej = __taylor_coeff(e, xk, alpha_k, j);
            if (ej.empty())
            {
                if (j < dk)
                {
                    xk_pow = xk_pow * xk_minus_alpha;
                    xk_pow.normalization();
                }
                continue;
            }

            // Step C: 递归多变量 Diophantine 求解 (GCL Algorithm 6.3)
            // eⱼ ∈ Z[x₁,...,x_{k-1}], 求 δᵢ 使得 eⱼ = Σ δᵢ · Ĝᵢ|_{xk=αk}
            auto deltas = __multivar_diophantine(
                ej, G_base, bezout_s, bezout_denom, v_factors,
                main_var, eval_vars, pa, 0);

            for (size_t i = 0; i < r; ++i)
            {
                if (deltas[i].empty()) continue;
                // Step D: Gᵢ += δᵢ · (xk - αk)^j
                auto update = deltas[i] * xk_pow;
                update.normalization();
                G[i] = G[i] + update;
                G[i].normalization();
            }
            prod_dirty = true;  // G 已更新, 下次迭代需重算 ∏Gᵢ

            // Step E (post-step): LC 校正
            for (size_t i = 0; i < r; ++i)
                __hensel_lc_correct(G[i], lc_tau[i], main_var);

            // 推进 (xk - αk)^j
            if (j < dk)
            {
                xk_pow = xk_pow * xk_minus_alpha;
                xk_pow.normalization();
            }
        }
    }

    // §6.2 + §6.7 多变量 Hensel 提升入口
    // 构造 Bézout 链, 初始化 Gᵢ, 逐变量提升
    template<class var_order>
    std::vector<polynomial_<ZZ, lex_<var_order>>>
    __multivar_hensel_lift(
        const polynomial_<ZZ, lex_<var_order>>& f_scaled,
        const std::vector<upolynomial_<ZZ>>& scaled_factors,
        const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_targets,
        const std::map<variable, ZZ>& eval_point,
        const variable& main_var)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        auto comp_ptr = f_scaled.comp_ptr();
        size_t r = scaled_factors.size();

        // ————————————————————————————————————————
        // §6.2 Bézout 链构造: Σ sᵢ·V̂ᵢ = denom
        // ————————————————————————————————————————
        std::vector<upolynomial_<ZZ>> bezout_s(r);
        ZZ denom(1);
        upolynomial_<ZZ> g_acc = scaled_factors[0];
        bezout_s[0] = upolynomial_<ZZ>({{umonomial(0), ZZ(1)}});  // s[0] = 1

        for (size_t i = 1; i < r; ++i)
        {
            upolynomial_<ZZ> alpha, beta;
            auto c_poly = polynomial_GCD(g_acc, scaled_factors[i], alpha, beta);
            // c_poly 应为常数 (pairwise coprime)
            assert(c_poly.size() == 1 && c_poly.front().first.deg() == 0);
            ZZ c_k = c_poly.front().second;

            // s[0..i-1] *= beta (无模约减)
            for (size_t j = 0; j < i; ++j)
            {
                bezout_s[j] = bezout_s[j] * beta;
                bezout_s[j].normalization();
            }
            // s[i] = denom * alpha
            auto denom_upoly = upolynomial_<ZZ>({{umonomial(0), denom}});
            bezout_s[i] = denom_upoly * alpha;
            bezout_s[i].normalization();

            denom *= c_k;
            g_acc = g_acc * scaled_factors[i];
            g_acc.normalization();
        }
        // 注: denom = ∏ c_k, 其中 c_k 是各步扩展 GCD 的内容
        // 对于整数多项式, denom 不一定为 ±1
        // 当 denom ≠ ±1 时, Diophantine 基本情形用模逆替代截断除法

        // ————————————————————————————————————————
        // §6.2.1 计算 p^a 模数 (仅当 denom ≠ ±1)
        // ————————————————————————————————————————
        ZZ pa(0);  // 哨兵: 0 = 精确除法 (denom == ±1)
        if (abs(denom) != ZZ(1))
        {
            // 系数界: ∏(||vi||₁ + 1) × |denom|
            ZZ B(1);
            for (auto& vi : scaled_factors)
            {
                ZZ vi_norm(0);
                for (auto& t : vi) vi_norm += abs(t.second);
                B *= vi_norm + ZZ(1);
            }
            B *= abs(denom);

            // 选素数 p 与 denom 和所有 lc(vi) 互素
            uint32_t p = 0;
            ZZ abs_denom = abs(denom);
            for (unsigned idx = 0; ; ++idx)
            {
                uint32_t candidate = boost::math::prime(idx);
                ZZ mod_d; ZZ::fdiv_r(mod_d, abs_denom, ZZ(candidate));
                if (mod_d == ZZ(0)) continue;
                bool bad = false;
                for (size_t i = 0; i < r; ++i)
                {
                    ZZ mod_lc; ZZ::fdiv_r(mod_lc, scaled_factors[i].front().second, ZZ(candidate));
                    if (mod_lc == ZZ(0)) { bad = true; break; }
                }
                if (bad) continue;
                p = candidate;
                break;
            }

            // pa = p^a > 2B
            pa = ZZ(1);
            ZZ two_B = ZZ(2) * B;
            while (pa <= two_B) pa *= ZZ(p);
        }

        // ————————————————————————————————————————
        // 初始化 Gᵢ: 将 vᵢ 转为 polynomial (仅含 main_var)
        // ————————————————————————————————————————
        std::vector<Poly> G(r, Poly(comp_ptr));
        for (size_t i = 0; i < r; ++i)
            poly_convert(scaled_factors[i], G[i], main_var);

        // δ = lc(vᵢ) (所有 vᵢ 的 lc 相同)
        ZZ delta = scaled_factors[0].front().second;

        // ————————————————————————————————————————
        // 确定提升变量顺序: x₂, x₃, ..., xₙ
        // ————————————————————————————————————————
        auto all_vars = get_variables(f_scaled);
        std::vector<variable> lift_vars;
        for (auto& [v, d] : all_vars)
            if (v != main_var)
                lift_vars.push_back(v);

        // ————————————————————————————————————————
        // 逐变量提升
        // ————————————————————————————————————————
        for (size_t k_idx = 0; k_idx < lift_vars.size(); ++k_idx)
        {
            variable xk = lift_vars[k_idx];
            ZZ alpha_k = eval_point.at(xk);

            // f_curr = f_scaled 代入 xk+1, ..., xn → α
            std::map<variable, ZZ> remaining_subs;
            for (size_t j = k_idx + 1; j < lift_vars.size(); ++j)
                remaining_subs[lift_vars[j]] = eval_point.at(lift_vars[j]);

            Poly f_curr = f_scaled;
            if (!remaining_subs.empty())
                f_curr = assign(f_curr, remaining_subs);

            int dk = degree(f_curr, xk);
            if (dk == 0) continue;

            // τᵢ 部分求值 (代入 xk+1, ..., xn)
            std::vector<Poly> lc_tau_curr(r, Poly(comp_ptr));
            for (size_t i = 0; i < r; ++i)
            {
                if (!remaining_subs.empty())
                    lc_tau_curr[i] = assign(lc_targets[i], remaining_subs);
                else
                    lc_tau_curr[i] = lc_targets[i];
            }

            // prev_eval: 前序已提升变量的求值点 {x₂→α₂,...,x_{k-1}→α_{k-1}}
            std::map<variable, ZZ> prev_eval;
            for (size_t j = 0; j < k_idx; ++j)
                prev_eval[lift_vars[j]] = eval_point.at(lift_vars[j]);

            __hensel_lift_one_var(
                f_curr, G, bezout_s, denom, pa, scaled_factors,
                delta, lc_tau_curr, main_var, xk, alpha_k, dk, prev_eval);
        }

        return G;
    }

    // ================================================================
    // §7. 多变量因式分解
    // ================================================================

    // §7.0 单项式 content 提取: 对每个变量计算所有项中的最小幂次，
    // 将公共变量幂次作为因子提取，返回除去公共幂次后的多项式。
    // 例：f = x²yz + x³y²z = x²yz(1 + xy)
    //   → var_factors = [(x, 2), (y, 1), (z, 1)], 返回 1 + xy
    template<class var_order>
    polynomial_<ZZ, lex_<var_order>>
    __extract_monomial_content(
        const polynomial_<ZZ, lex_<var_order>>& f,
        std::vector<std::pair<variable, int64_t>>& var_factors)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        var_factors.clear();

        if (f.empty()) return f;

        // Step 1: 各变量取所有项中的最小幂次
        std::map<variable, int64_t> min_deg;
        bool first_term = true;

        for (const auto& [mono, coeff] : f)
        {
            if (first_term)
            {
                for (const auto& [var, deg] : mono)
                    min_deg[var] = deg;
                first_term = false;
            }
            else
            {
                std::set<variable> present;
                for (const auto& [var, deg] : mono)
                {
                    present.insert(var);
                    auto it = min_deg.find(var);
                    if (it != min_deg.end())
                        it->second = std::min(it->second, deg);
                }
                // 不在本项中的变量 → min 降为 0 → 移除
                auto it = min_deg.begin();
                while (it != min_deg.end())
                {
                    if (present.find(it->first) == present.end())
                        it = min_deg.erase(it);
                    else
                        ++it;
                }
            }
        }

        // Step 2: 移除 min_deg 为 0 的变量
        for (auto it = min_deg.begin(); it != min_deg.end(); )
        {
            if (it->second == 0)
                it = min_deg.erase(it);
            else
                ++it;
        }

        if (min_deg.empty()) return f;  // 无公共变量幂

        // Step 3: 构造提取后的多项式（每项的单项式减去 min_deg）
        Poly result(f.comp_ptr());
        for (const auto& [mono, coeff] : f)
        {
            typename Poly::monomial_type new_mono;
            for (const auto& [var, deg] : mono)
            {
                auto it = min_deg.find(var);
                int64_t subtract = (it != min_deg.end()) ? it->second : 0;
                if (deg > subtract)
                    new_mono.push_back({var, deg - subtract});
            }
            result.data().push_back({std::move(new_mono), coeff});
        }
        result.normalization();

        // Step 4: 记录提取的变量因子
        for (const auto& [var, deg] : min_deg)
            var_factors.push_back({var, deg});

        return result;
    }

    // §7.1 Wang 核心: 对无平方本原多变量多项式执行 Wang 算法
    template<class var_order>
    std::vector<std::pair<polynomial_<ZZ, lex_<var_order>>, uint64_t>>
    __wang_core(const polynomial_<ZZ, lex_<var_order>>& g)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        auto comp_ptr = g.comp_ptr();

        auto all_vars_list = get_variables(g);
        std::vector<std::pair<variable, int64_t>> all_vars(
            all_vars_list.begin(), all_vars_list.end());

        // 交错轮换主变量: 每个主变量分配 BATCH_SIZE 个求值点,
        // 用完后换下一个主变量, 一轮结束后扩展配额继续.
        // 避免在某个"坏"主变量上浪费过多时间.
        // 循环直到所有主变量都被标记为 dead (不可约启发式).
        const int BATCH_SIZE = 200;

        // 筛选可用主变量 (deg > 1)
        std::vector<variable> main_vars;
        for (auto& [v, d] : all_vars)
            if (degree(g, v) > 1)
                main_vars.push_back(v);

        if (main_vars.empty())
            return {{g, 1}};

        // 每个主变量的状态
        std::vector<int> var_skip(main_vars.size(), 0);         // 当前 skip 位置
        std::vector<bool> var_dead(main_vars.size(), false);    // 已证明不可约

        // 辅助 lambda: 规范化原始部分 (pp + 正首项)
        auto normalize_factor = [](Poly& h)
        {
            h = pp(h);
            if (!h.empty() && h.front().second < 0)
                for (auto& term : h.data())
                    term.second = -term.second;
        };

        // 终止性: 可约多项式至少存在一个好求值点使 Hensel 提升成功 (→ return verified);
        // 不可约多项式的所有主变量最终被标记为 dead (→ break via all_dead).
        for (;;)
        {
            for (size_t vi = 0; vi < main_vars.size(); ++vi)
            {
                if (var_dead[vi]) continue;
                variable x1 = main_vars[vi];
                int batch_end = var_skip[vi] + BATCH_SIZE;

                for (int skip = var_skip[vi]; skip < batch_end; ++skip)
                {
                    auto eval = __select_eval_point(g, x1, skip);
                    auto f0 = assign(g, eval);
                    upolynomial_<ZZ> f0_upoly;
                    poly_convert(f0, f0_upoly);
                    auto uni_fac = factorize(f0_upoly);

                    std::vector<upolynomial_<ZZ>> uni_factors;
                    for (auto& [fi, ei] : uni_fac.factors)
                        uni_factors.push_back(fi);

                    if (uni_factors.size() <= 1)
                    {
                        // f(x₁, α) 在满足条件 (a)(b) 的 α 处不可约
                        // → 数学上证明 g 关于 x₁ 无正次因子分解
                        // 标记此主变量已证明不可约, 换下一个主变量
                        var_dead[vi] = true;
                        break;
                    }

                    auto lc_result = __wang_leading_coeff(
                        g, uni_factors, eval, x1, uni_fac.content);
                    if (!lc_result.success)
                        continue;

                    auto mv_factors = __multivar_hensel_lift(
                        lc_result.f_scaled, lc_result.scaled_factors,
                        lc_result.lc_targets, eval, x1);

                    // 因子重组 (Zassenhaus 子集枚举, 同 §7.5 单变量版本)
                    // 从 s=1 递增枚举 Hensel 因子子集，找最小整除子集。
                    // 不可约性保证: 若返回因子 F 可约 (F=A·B)，则 A,B 各对应
                    // Hensel 因子的真子集 (模求值回到不可约单变量因子，唯一分解)，
                    // 即存在更小的 s 使子集乘积整除 g，与 F 在 s 处被找到矛盾。
                    // 预处理: 规范化所有 Hensel 因子
                    std::vector<Poly> normed(mv_factors.size(), Poly(comp_ptr));
                    std::vector<size_t> active_idx;
                    for (size_t fi = 0; fi < mv_factors.size(); ++fi)
                    {
                        normed[fi] = mv_factors[fi];
                        normalize_factor(normed[fi]);
                        if (!normed[fi].empty() && !is_number(normed[fi]))
                            active_idx.push_back(fi);
                    }

                    std::vector<std::pair<Poly, uint64_t>> verified;
                    Poly g_remaining = g;

                    // T: 可用因子在 active_idx 中的下标
                    std::vector<int> mv_T(active_idx.size());
                    std::iota(mv_T.begin(), mv_T.end(), 0);

                    auto mv_next_combination = [](std::vector<int>& idx, int n) -> bool {
                        int sz = (int)idx.size();
                        int i = sz - 1;
                        while (i >= 0 && idx[i] == n - sz + i)
                            --i;
                        if (i < 0) return false;
                        ++idx[i];
                        for (int j = i + 1; j < sz; ++j)
                            idx[j] = idx[j - 1] + 1;
                        return true;
                    };

                    if (!mv_T.empty())
                    {
                        bool found;
                        int s = 1;
                        while (2 * s <= (int)mv_T.size())
                        {
                            found = false;
                            int n_cur = (int)mv_T.size();
                            if (s > n_cur / 2) break;

                            std::vector<int> idx(s);
                            std::iota(idx.begin(), idx.end(), 0);

                            do {
                                // 计算子集乘积
                                Poly prod = normed[active_idx[mv_T[idx[0]]]];
                                for (int k = 1; k < s; ++k)
                                {
                                    prod = prod * normed[active_idx[mv_T[idx[k]]]];
                                    prod.normalization();
                                }
                                normalize_factor(prod);

                                if (!prod.empty() && !is_number(prod))
                                {
                                    Poly q(comp_ptr), rem(comp_ptr);
                                    pair_vec_div(q.data(), rem.data(),
                                                 g_remaining.data(), prod.data(), g.comp());
                                    if (rem.empty() && !q.empty())
                                    {
                                        verified.push_back({std::move(prod), 1});
                                        g_remaining = std::move(q);
                                        normalize_factor(g_remaining);
                                        for (int j = s - 1; j >= 0; --j)
                                            mv_T.erase(mv_T.begin() + idx[j]);
                                        found = true;
                                        break;
                                    }
                                }
                            } while (mv_next_combination(idx, n_cur));

                            if (found)
                            {
                                s = 1;
                                continue;
                            }
                            ++s;
                        }
                    }

                    // 剩余部分作为最后一个因子
                    if (!g_remaining.empty() && !is_number(g_remaining))
                    {
                        auto h = g_remaining;
                        normalize_factor(h);
                        if (!is_number(h))
                            verified.push_back({std::move(h), 1});
                    }

                    if (verified.size() >= 2)
                        return verified;
                }
                var_skip[vi] = batch_end;
            }

            // 检查是否所有主变量都已 dead
            bool all_dead = true;
            for (size_t vi = 0; vi < main_vars.size(); ++vi)
                if (!var_dead[vi]) { all_dead = false; break; }
            if (all_dead) break;
        }

        return {{g, 1}};
    }

    // §7.1 多变量因式分解入口
    template<class var_order>
    factorization<polynomial_<ZZ, lex_<var_order>>>
    __factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f_input)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        factorization<Poly> result;
        result.content = ZZ(1);

        auto sqf = squarefreefactorize(f_input);

        for (auto& [gk, mk] : sqf)
        {
            if (is_number(gk))
            {
                ZZ c = gk.front().second;
                for (uint64_t e = 0; e < mk; ++e)
                    result.content *= c;
                continue;
            }

            auto gk_vars = get_variables(gk);
            if (gk_vars.size() <= 1)
            {
                auto sub = factorize(gk);
                for (uint64_t e = 0; e < mk; ++e)
                    result.content *= sub.content;
                for (auto& [fi, ei] : sub.factors)
                    result.factors.push_back({fi, ei * mk});
            }
            else
            {
                // 预处理：提取单项式 content（公共变量幂次）
                std::vector<std::pair<variable, int64_t>> mono_factors;
                Poly gk_reduced = __extract_monomial_content(gk, mono_factors);

                // 将提取的变量幂次作为因子加入结果
                for (auto& [var, vdeg] : mono_factors)
                {
                    Poly var_poly(gk.comp_ptr());
                    typename Poly::monomial_type mono;
                    mono.push_back({var, 1});
                    var_poly.data().push_back({std::move(mono), ZZ(1)});
                    var_poly.normalization();
                    result.factors.push_back({std::move(var_poly), static_cast<uint64_t>(vdeg) * mk});
                }

                // 检查提取后是否还需要多变量分解
                auto reduced_vars = get_variables(gk_reduced);
                if (reduced_vars.size() <= 1)
                {
                    // 降为单变量或常数
                    auto sub = factorize(gk_reduced);
                    for (uint64_t e = 0; e < mk; ++e)
                        result.content *= sub.content;
                    for (auto& [fi, ei] : sub.factors)
                        result.factors.push_back({fi, ei * mk});
                }
                else
                {
                    // 确保传给 __wang_core 的多项式 lc > 0
                    Poly gk_pos = gk_reduced;
                    bool negated = false;
                    if (!gk_pos.empty() && gk_pos.front().second < 0)
                    {
                        for (auto& term : gk_pos.data())
                            term.second = -term.second;
                        negated = true;
                    }
                    auto wang_factors = __wang_core(gk_pos);
                    for (auto& [fi, ei] : wang_factors)
                        result.factors.push_back({fi, ei * mk});
                    if (negated)
                    {
                        // (-1)^mk 吸收到 content
                        if (mk % 2 == 1)
                            result.content = -result.content;
                    }
                }
            }
        }

        // 确保所有因子 lc > 0
        for (auto& [fac, mult] : result.factors)
        {
            if (!fac.empty() && fac.front().second < 0)
            {
                for (auto& term : fac.data())
                    term.second = -term.second;
                if (mult % 2 != 0)
                    result.content = -result.content;
            }
        }

        std::sort(result.factors.begin(), result.factors.end(),
            [](const auto& a, const auto& b) {
                return degree(a.first) < degree(b.first);
            });

#ifndef NDEBUG
        {
            Poly check(f_input.comp_ptr());
            check.data().push_back({{}, result.content});
            check.normalization();
            for (const auto& [fi, ei] : result.factors)
                for (uint64_t e = 0; e < ei; ++e)
                {
                    check = check * fi;
                    check.normalization();
                }
            assert(check == f_input);
        }
#endif

        return result;
    }

} // namespace clpoly

#endif // CLPOLY_POLYNOMIAL_FACTORIZE_WANG_HH
