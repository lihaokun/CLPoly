/**
 * @file polynomial_factorize.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 多项式因式分解
 *
 */
#ifndef CLPOLY_POLYNOMIAL_FACTORIZE_HH
#define CLPOLY_POLYNOMIAL_FACTORIZE_HH

#include <clpoly/polynomial.hh>
#include <clpoly/polynomial_gcd.hh>
#include <clpoly/upolynomial.hh>
#include <boost/math/special_functions/prime.hpp>
#include <random>
#include <algorithm>
#include <stdexcept>

namespace clpoly{

    // ================================================================
    // §8.1 factorization 返回类型
    // ================================================================
    template<class Poly>
    struct factorization {
        typename Poly::coeff_type content;
        std::vector<std::pair<Poly, uint64_t>> factors;
    };

    // Helper: 构造 Zp(int64_t, uint32_t) 避免重载歧义
    inline Zp __make_zp(int64_t val, uint32_t p) { return Zp(val, p); }

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

    // §4.1.5-6 polynomial_GCD (普通+扩展) 已迁移到 polynomial_gcd.hh

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
        const variable& main_var)
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
                auto f0 = assign(f, alpha);
                if (!is_squarefree(f0)) continue;

                // 条件 (d): lc 因子求值两两互素
                if (!lc_irr_factors.empty())
                {
                    std::vector<ZZ> vals;
                    bool skip = false;
                    for (auto& lj : lc_irr_factors)
                    {
                        auto lj_eval = assign(lj, alpha);
                        if (lj_eval.empty()) { skip = true; break; }
                        ZZ v = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
                        if (v == 0) { skip = true; break; }
                        vals.push_back(abs(v));
                    }
                    if (skip) continue;
                    // 条件 (d'): |lⱼ(α)| ≥ 2, 避免退化赋值
                    for (auto& v : vals)
                        if (v <= ZZ(1)) { skip = true; break; }
                    if (skip) continue;
                    // 两两互素
                    bool pairwise_coprime = true;
                    for (size_t i = 0; i < vals.size() && pairwise_coprime; ++i)
                        for (size_t j = i + 1; j < vals.size() && pairwise_coprime; ++j)
                            if (gcd(vals[i], vals[j]) != ZZ(1))
                                pairwise_coprime = false;
                    if (!pairwise_coprime) continue;
                }

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
        const variable& main_var)
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

            // wᵢ = lc(uᵢ) — 数值跟踪
            std::vector<ZZ> w(r);
            for (size_t i = 0; i < r; ++i)
                w[i] = univar_factors[i].front().second;

            for (auto& [lj, ej] : lc_factors)
            {
                // zⱼ = |lⱼ(α)|^eⱼ
                auto lj_eval = assign(lj, eval_point);
                ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
                if (lj_val == 0) return result;  // 不应发生（条件 d 保证）
                ZZ zj(1);
                for (uint64_t e = 0; e < ej; ++e)
                    zj *= lj_val;

                // 找唯一的 uᵢ 使得 zⱼ | wᵢ
                std::vector<size_t> candidates;
                for (size_t i = 0; i < r; ++i)
                    if (w[i] % zj == ZZ(0))
                        candidates.push_back(i);

                if (candidates.size() != 1) return result;  // FAIL

                size_t idx = candidates[0];
                // σᵢ *= lⱼ^eⱼ
                Poly lj_power = lj;
                for (uint64_t e = 1; e < ej; ++e)
                {
                    lj_power = lj_power * lj;
                    lj_power.normalization();
                }
                sigma[idx] = sigma[idx] * lj_power;
                sigma[idx].normalization();
                w[idx] /= zj;
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

    // §6.5 辅助: 多项式每个系数除以整数 d (精确整除)
    template<class var_order>
    polynomial_<ZZ, lex_<var_order>> __poly_exact_div_zz(
        const polynomial_<ZZ, lex_<var_order>>& f,
        const ZZ& d)
    {
        if (d == ZZ(1)) return f;
        polynomial_<ZZ, lex_<var_order>> result(f.comp_ptr());
        result.reserve(f.size());
        for (auto& term : f)
        {
            ZZ q = term.second / d;
            if (q != 0)
                result.push_back({term.first, q});
        }
        return result;
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

        // (xk - αk)^j, 从 j=1 开始
        Poly xk_pow = xk_minus_alpha;  // (xk - αk)^1

        for (int j = 1; j <= dk; ++j)
        {
            // Step A: e = f_curr - ∏Gᵢ
            Poly prod = G[0];
            for (size_t i = 1; i < r; ++i)
            {
                prod = prod * G[i];
                prod.normalization();
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

            // Step C: MDP — 先将 eⱼ 求值到单变量，再用 Bézout 求解
            // eⱼ ∈ Z[x₁,...,x_{k-1}], 求值前序变量 → ē_j ∈ Z[x₁]
            Poly ej_eval = ej;
            if (!prev_eval.empty())
                ej_eval = assign(ej, prev_eval);

            // 转为 upolynomial 进行单变量 MDP
            upolynomial_<ZZ> ej_upoly;
            poly_convert(ej_eval, ej_upoly);

            for (size_t i = 0; i < r; ++i)
            {
                auto si_ej = bezout_s[i] * ej_upoly;
                si_ej.normalization();

                int64_t deg_si_ej = get_deg(si_ej);
                int64_t deg_vi = get_deg(v_factors[i]);
                int64_t ki = std::max(deg_si_ej - deg_vi + 1, (int64_t)0);

                // prem 在 upolynomial 上: lc(vi)^ki · si_ej = Q·vi + prem
                // 因为 upolynomial 没有 prem, 用 pair_vec_div 乘以 lc^ki 后除
                ZZ lc_vi = v_factors[i].front().second;  // = δ
                ZZ lc_pow(1);
                for (int64_t e = 0; e < ki; ++e)
                    lc_pow *= lc_vi;

                // scaled = lc^ki · si_ej
                upolynomial_<ZZ> scaled;
                if (lc_pow != ZZ(1))
                {
                    upolynomial_<ZZ> lc_poly({{umonomial(0), lc_pow}});
                    scaled = lc_poly * si_ej;
                    scaled.normalization();
                }
                else
                    scaled = si_ej;

                // pair_vec_div: scaled = Q · vi + rem
                upolynomial_<ZZ> q_upoly, rem_upoly;
                pair_vec_div(q_upoly.data(), rem_upoly.data(),
                             scaled.data(), v_factors[i].data(), v_factors[i].comp());

                // δ̄ᵢ = rem / (δ^kᵢ · denom)
                ZZ divisor = bezout_denom * lc_pow;
                upolynomial_<ZZ> delta_i_upoly;
                delta_i_upoly.reserve(rem_upoly.size());
                for (auto& term : rem_upoly)
                {
                    ZZ coeff = term.second / divisor;
                    if (coeff != 0)
                        delta_i_upoly.push_back({term.first, coeff});
                }

                // 转为 multivariate polynomial
                Poly delta_i(comp_ptr);
                poly_convert(delta_i_upoly, delta_i, main_var);

                // Step D: Gᵢ += δ̄ᵢ · (xk - αk)^j
                auto update = delta_i * xk_pow;
                update.normalization();
                G[i] = G[i] + update;
                G[i].normalization();
            }

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
                f_curr, G, bezout_s, denom, scaled_factors,
                delta, lc_tau_curr, main_var, xk, alpha_k, dk, prev_eval);
        }

        return G;
    }

    // §7 前向声明
    template<class var_order>
    factorization<polynomial_<ZZ, lex_<var_order>>>
    __factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f_input);

    // §4.1.7 模幂: base^exp mod modpoly in Zp[x]
    inline upolynomial_<Zp> __upoly_powmod(
        const upolynomial_<Zp>& base,
        const ZZ& exp,
        const upolynomial_<Zp>& modpoly)
    {
        assert(!modpoly.empty());
        uint32_t p = modpoly.front().second.prime();
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
        uint32_t p,
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
        uint32_t p = f.front().second.prime();
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
        uint32_t p = f.front().second.prime();

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
    inline upolynomial_<Zp> __upoly_subtract_x(const upolynomial_<Zp>& h, uint32_t p)
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
    inline upolynomial_<Zp> __upoly_subtract_one(const upolynomial_<Zp>& h, uint32_t p)
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
        uint32_t p = f.front().second.prime();
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

        uint32_t p = f.front().second.prime();
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
        uint32_t p = f.front().second.prime();

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
        assert(r <= 64);  // 位掩码限制

        // T: 可用因子位掩码 (第 i 位 = 1 表示 lifted[i] 可用)
        uint64_t T = (r == 64) ? ~0ULL : ((1ULL << r) - 1);
        upolynomial_<ZZ> f_star = f;
        std::vector<upolynomial_<ZZ>> result;

        // popcount helper
        auto popcount = [](uint64_t x) -> int {
            int c = 0;
            while (x) { c += x & 1; x >>= 1; }
            return c;
        };

        // 从位掩码提取子集下标
        auto mask_to_subset = [](uint64_t mask) -> std::vector<size_t> {
            std::vector<size_t> sub;
            for (size_t i = 0; i < 64; ++i)
                if (mask & (1ULL << i)) sub.push_back(i);
            return sub;
        };

        // T 中的因子数
        auto T_count = [&popcount](uint64_t mask) -> int {
            return popcount(mask);
        };

        bool found;
        int s = 1;
        while (2 * s <= T_count(T))
        {
            found = false;

            // Gosper's hack: 枚举 T 中大小为 s 的子集
            // 先枚举 r 位中恰好 s 个 1 的所有掩码, 再过滤只含 T 中的位
            // 更简单: 枚举 T 的子集中大小为 s 的

            // 收集 T 中的活跃位
            std::vector<size_t> active;
            for (size_t i = 0; i < r; ++i)
                if (T & (1ULL << i)) active.push_back(i);
            int n_active = (int)active.size();

            if (s > n_active / 2) break;

            // 枚举 active 的大小为 s 的子集 (用 Gosper's hack on n_active 位)
            uint64_t sub_mask = (1ULL << s) - 1;
            uint64_t sub_limit = (1ULL << n_active);

            while (sub_mask < sub_limit)
            {
                // 将 sub_mask 映射回原始因子索引
                std::vector<size_t> S_idx;
                uint64_t S_bits = 0;
                for (int j = 0; j < n_active; ++j)
                {
                    if (sub_mask & (1ULL << j))
                    {
                        S_idx.push_back(active[j]);
                        S_bits |= (1ULL << active[j]);
                    }
                }

                ZZ lc_fstar = f_star.front().second;

                // === 剪枝 1: 首项系数检查 ===
                ZZ lc_prod = lc_fstar;
                for (size_t idx : S_idx)
                    lc_prod *= lifted[idx].front().second;
                lc_prod = __symmetric_mod(lc_prod, m);
                ZZ lc_sq = lc_fstar * lc_fstar;
                // lc_prod 是否整除 lc(f*)²?
                if (lc_prod != ZZ(0))
                {
                    ZZ rem;
                    ZZ::fdiv_r(rem, lc_sq, lc_prod);
                    if (rem != ZZ(0))
                        goto next_subset;
                }

                {
                    // === 剪枝 2: 常数项检查 ===
                    ZZ c_prod = lc_fstar;
                    for (size_t idx : S_idx)
                        c_prod *= __upoly_const_term(lifted[idx]);
                    c_prod = __symmetric_mod(c_prod, m);
                    ZZ fstar_const = lc_fstar * __upoly_const_term(f_star);
                    if (c_prod != ZZ(0))
                    {
                        ZZ rem2;
                        ZZ::fdiv_r(rem2, fstar_const, c_prod);
                        if (rem2 != ZZ(0))
                            goto next_subset;
                    }

                    {
                        // === 完整验证: 试除 ===
                        auto g = __subset_product_mod(lifted, S_idx, lc_fstar, m);

                        // lc(f*) * f*
                        upolynomial_<ZZ> lc_f_star_poly = f_star;
                        for (auto& term : lc_f_star_poly)
                            term.second *= lc_fstar;
                        lc_f_star_poly.normalization();

                        // 试除: g | lc(f*)*f* in Z[x]?
                        upolynomial_<ZZ> q_trial, r_trial;
                        pair_vec_div(q_trial.data(), r_trial.data(),
                                     lc_f_star_poly.data(), g.data(), f_star.comp());

                        if (r_trial.empty())
                        {
                            // 找到真因子
                            auto [c_g, pp_g] = __upoly_primitive(std::move(g));
                            result.push_back(std::move(pp_g));

                            auto [c_q, pp_q] = __upoly_primitive(std::move(q_trial));
                            f_star = std::move(pp_q);

                            uint64_t T_comp_bits = T & ~S_bits;
                            T = T_comp_bits;
                            found = true;
                            break;
                        }
                    }
                }

                next_subset:
                // Gosper's hack: 下一个同样 popcount 的掩码
                if (sub_mask == 0) break;
                uint64_t c_ = sub_mask & (-sub_mask);
                uint64_t r_ = sub_mask + c_;
                sub_mask = (((r_ ^ sub_mask) >> 2) / c_) | r_;
            }

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
        size_t max_tries = 5;
        std::mt19937 rng(42);

        for (size_t idx = 0, tried = 0; tried < max_tries; ++idx)
        {
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

            // 因子数 > deg/2 时扩展尝试到 20 个
            if (best_count > (size_t)(deg_f / 2) && tried == 5)
                max_tries = 20;
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

        return result;
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

        return result;
    }

    // ================================================================
    // §7. 多变量因式分解
    // ================================================================

    // §7.1 Wang 核心: 对无平方本原多变量多项式执行 Wang 算法
    template<class var_order>
    std::vector<std::pair<polynomial_<ZZ, lex_<var_order>>, uint64_t>>
    __wang_core(const polynomial_<ZZ, lex_<var_order>>& g)
    {
        using Poly = polynomial_<ZZ, lex_<var_order>>;
        auto comp_ptr = g.comp_ptr();

        auto all_vars = get_variables(g);
        variable x1 = all_vars.front().first;

        const int MAX_RETRY = 10;

        for (int retry = 0; retry < MAX_RETRY; ++retry)
        {
            auto eval = __select_eval_point(g, x1);

            auto f0 = assign(g, eval);
            upolynomial_<ZZ> f0_upoly;
            poly_convert(f0, f0_upoly);
            auto uni_fac = factorize(f0_upoly);

            std::vector<upolynomial_<ZZ>> uni_factors;
            for (auto& [fi, ei] : uni_fac.factors)
                uni_factors.push_back(fi);

            if (uni_factors.size() <= 1)
                return {{g, 1}};

            auto lc_result = __wang_leading_coeff(g, uni_factors, eval, x1);
            if (!lc_result.success)
                continue;

            auto mv_factors = __multivar_hensel_lift(
                lc_result.f_scaled, lc_result.scaled_factors,
                lc_result.lc_targets, eval, x1);

            // 试除验证 + 去缩放
            std::vector<std::pair<Poly, uint64_t>> verified;
            Poly g_remaining = g;

            for (auto& Gi : mv_factors)
            {
                auto h = pp(Gi);
                if (!h.empty() && h.front().second < 0)
                {
                    for (auto& term : h.data())
                        term.second = -term.second;
                }

                Poly q(comp_ptr), rem(comp_ptr);
                pair_vec_div(q.data(), rem.data(),
                             g_remaining.data(), h.data(), g.comp());
                if (rem.empty() && !q.empty())
                {
                    verified.push_back({std::move(h), 1});
                    g_remaining = std::move(q);
                }
            }

            if (!g_remaining.empty() && !is_number(g_remaining))
            {
                auto h = pp(g_remaining);
                if (!h.empty() && h.front().second < 0)
                {
                    for (auto& term : h.data())
                        term.second = -term.second;
                }
                if (!is_number(h))
                    verified.push_back({std::move(h), 1});
            }

            if (!verified.empty())
                return verified;
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
                // 确保传给 __wang_core 的多项式 lc > 0
                Poly gk_pos = gk;
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

        // 确保所有因子 lc > 0
        for (auto& [fac, mult] : result.factors)
        {
            if (!fac.empty() && fac.front().second < 0)
            {
                for (auto& term : fac.data())
                    term.second = -term.second;
                result.content = -result.content;
            }
        }

        std::sort(result.factors.begin(), result.factors.end(),
            [](const auto& a, const auto& b) {
                return degree(a.first) < degree(b.first);
            });

        return result;
    }

} // namespace clpoly
#endif
