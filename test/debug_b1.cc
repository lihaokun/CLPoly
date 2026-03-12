/**
 * @file debug_b1.cc
 * @brief Diagnostic test for B1 factorization bug
 *
 * f = f1*f2*f3*f4 where:
 *   f1 = y^2 - y - 2 = (y-2)(y+1)
 *   f2 = -2x^2 + xy - y
 *   f3 = -2x^2 - 2y^2 + 3
 *   f4 = x^2 - 3y^2 + y
 *
 * Expected: 5 irreducible factors (y-2)(y+1)(f2)(f3)(f4)
 * Bug: returns 4 factors, with f2*f4 fused.
 *
 * Build:
 *   g++ -g -D DEBUG -I/home/haokun/projects/CLPoly /tmp/debug_b1.cc \
 *       -o /tmp/debug_b1 /home/haokun/projects/CLPoly/lib/debug/libclpoly.a \
 *       -lgmpxx -lgmp
 */

#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>
#include <iomanip>

using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;

// Helper to print a polynomial nicely
void print_poly(const char* label, const PolyZZ& p) {
    std::cout << label << " = " << p << std::endl;
}

void print_upoly(const char* label, const upolynomial_<ZZ>& p) {
    std::cout << label << " = ";
    for (size_t i = 0; i < p.size(); ++i) {
        if (i > 0) std::cout << " + ";
        std::cout << p[i].second << "*x^" << p[i].first.deg();
    }
    std::cout << std::endl;
}

int main()
{
    variable x("x"), y("y");

    // ============================================================
    // Step 0: Construct f = f1 * f2 * f3 * f4
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Step 0: Construct f = f1*f2*f3*f4\n";
    std::cout << "========================================\n";

    PolyZZ f1_lex, f2_lex, f3_lex, f4_lex;
    {
        polynomial_ZZ f1 = pow(y,2) - y - ZZ(2);
        polynomial_ZZ f2 = ZZ(-2)*pow(x,2) + x*y - y;
        polynomial_ZZ f3 = ZZ(-2)*pow(x,2) - ZZ(2)*pow(y,2) + ZZ(3);
        polynomial_ZZ f4 = pow(x,2) - ZZ(3)*pow(y,2) + y;

        poly_convert(f1, f1_lex);
        poly_convert(f2, f2_lex);
        poly_convert(f3, f3_lex);
        poly_convert(f4, f4_lex);
    }

    print_poly("f1", f1_lex);
    print_poly("f2", f2_lex);
    print_poly("f3", f3_lex);
    print_poly("f4", f4_lex);

    PolyZZ f_lex = f1_lex * f2_lex;
    f_lex.normalization();
    f_lex = f_lex * f3_lex;
    f_lex.normalization();
    f_lex = f_lex * f4_lex;
    f_lex.normalization();

    print_poly("f = f1*f2*f3*f4", f_lex);
    std::cout << std::endl;

    // ============================================================
    // Step 1: squarefreefactorize(f)
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Step 1: squarefreefactorize(f)\n";
    std::cout << "========================================\n";

    auto sqf = squarefreefactorize(f_lex);
    std::cout << "Number of squarefree factors: " << sqf.size() << "\n";
    for (size_t i = 0; i < sqf.size(); ++i) {
        std::cout << "  sqf[" << i << "] (mult=" << sqf[i].second << "): "
                  << sqf[i].first << "\n";
        auto vars = get_variables(sqf[i].first);
        std::cout << "    variables: ";
        for (auto& [v, d] : vars)
            std::cout << v << "(deg=" << d << ") ";
        if (is_number(sqf[i].first))
            std::cout << "[constant]";
        else if (vars.size() == 1)
            std::cout << "[univariate]";
        else
            std::cout << "[multivariate]";
        std::cout << "\n";
    }
    std::cout << std::endl;

    // ============================================================
    // Step 2: For each multivariate sqf factor, simulate __wang_core
    // ============================================================
    std::cout << "========================================\n";
    std::cout << "Step 2: Simulate __wang_core for multivariate factors\n";
    std::cout << "========================================\n";

    for (size_t si = 0; si < sqf.size(); ++si) {
        auto& g = sqf[si].first;
        auto vars = get_variables(g);
        if (is_number(g) || vars.size() <= 1) continue;

        std::cout << "\n--- Processing sqf factor " << si << ": " << g << "\n";

        // Extract monomial content
        std::vector<std::pair<variable, int64_t>> mono_factors;
        // Replicate __extract_monomial_content inline
        {
            // Find minimum degree of each variable across all terms
            std::map<variable, int64_t> min_degs;
            bool first = true;
            for (const auto& term : g) {
                std::map<variable, int64_t> term_degs;
                for (const auto& vp : term.first)
                    term_degs[vp.first] = vp.second;
                if (first) {
                    min_degs = term_degs;
                    first = false;
                } else {
                    // For each var in min_degs, if not in this term, min = 0
                    for (auto it = min_degs.begin(); it != min_degs.end(); ) {
                        auto jt = term_degs.find(it->first);
                        if (jt == term_degs.end())
                            it = min_degs.erase(it);
                        else {
                            it->second = std::min(it->second, jt->second);
                            ++it;
                        }
                    }
                }
            }
            for (auto& [v, d] : min_degs)
                if (d > 0)
                    mono_factors.push_back({v, d});
        }

        if (!mono_factors.empty()) {
            std::cout << "  Monomial content extracted: ";
            for (auto& [v, d] : mono_factors)
                std::cout << v << "^" << d << " ";
            std::cout << "\n";
        } else {
            std::cout << "  No monomial content.\n";
        }

        // Ensure lc > 0
        PolyZZ g_pos = g;
        if (!g_pos.empty() && g_pos.front().second < 0) {
            for (auto& term : g_pos.data())
                term.second = -term.second;
            std::cout << "  Negated to make lc > 0: " << g_pos << "\n";
        }

        // For each candidate main variable
        auto all_vars = get_variables(g_pos);
        for (auto& [main_v, main_d] : all_vars) {
            if (main_d <= 1) continue;

            std::cout << "\n  === Main variable: " << main_v
                      << " (deg=" << main_d << ") ===\n";

            // leadcoeff(g_pos, main_v)
            auto L = leadcoeff(g_pos, main_v);
            std::cout << "  lc(g, " << main_v << ") = " << L << "\n";

            // If lc is non-constant, factorize it
            if (!is_number(L)) {
                auto lc_fac = factorize(L);
                std::cout << "  lc factorization: content=" << lc_fac.content << ", factors:\n";
                for (auto& [lf, le] : lc_fac.factors)
                    std::cout << "    (" << lf << ")^" << le << "\n";
            }

            // Try first few eval points
            int max_points = 5;
            uint64_t mtshl_p = UINT64_C(18446744073709551557);  // 2^64 - 59

            // Check if mtshl_p divides lc
            if (!is_number(L)) {
                auto all_div = [&](uint64_t p) {
                    for (const auto& term : L)
                        if (term.second.fdiv_ui(p) != 0) return false;
                    return true;
                };
                while (all_div(mtshl_p))
                    mtshl_p = prev_prime_64(mtshl_p);
            }

            // Also factorize lc for condition (d) check
            std::vector<PolyZZ> lc_irr_factors;
            if (!is_number(L)) {
                auto lc_fac = factorize(L);
                for (auto& [lj, ej] : lc_fac.factors)
                    lc_irr_factors.push_back(lj);
            }

            int found_points = 0;
            for (int skip = 0; skip < 50 && found_points < max_points; ++skip) {
                // Use __select_eval_point via replicated logic
                auto eval = __select_eval_point(g_pos, main_v, skip, mtshl_p);

                std::cout << "\n  --- Eval point #" << found_points << " (skip=" << skip << "): ";
                for (auto& [v, val] : eval)
                    std::cout << v << "=" << val << " ";
                std::cout << "\n";

                // Evaluate: f0 = g_pos(main_v, alpha)
                auto f0 = assign(g_pos, eval);
                std::cout << "  f0 = g(x1, alpha) = " << f0 << "\n";

                // Convert to upoly and factor
                upolynomial_<ZZ> f0_upoly;
                poly_convert(f0, f0_upoly);
                print_upoly("  f0 as upoly", f0_upoly);

                auto uni_fac = factorize(f0_upoly);
                std::cout << "  Univariate factorization:\n";
                std::cout << "    content = " << uni_fac.content << "\n";
                std::cout << "    # factors = " << uni_fac.factors.size() << "\n";
                for (size_t fi = 0; fi < uni_fac.factors.size(); ++fi) {
                    std::cout << "    factor[" << fi << "] (mult="
                              << uni_fac.factors[fi].second << "): ";
                    auto& uf = uni_fac.factors[fi].first;
                    for (size_t ti = 0; ti < uf.size(); ++ti) {
                        if (ti > 0) std::cout << " + ";
                        std::cout << uf[ti].second << "*x^" << uf[ti].first.deg();
                    }
                    std::cout << "\n";
                }

                // Check if LC correction would succeed
                std::vector<upolynomial_<ZZ>> uni_factors;
                for (auto& [fi, ei] : uni_fac.factors)
                    uni_factors.push_back(fi);

                if (uni_factors.size() <= 1) {
                    std::cout << "  -> IRREDUCIBLE image (1 or 0 factors)\n";
                    std::cout << "  -> __wang_core would mark " << main_v
                              << " as DEAD for this main variable\n";
                } else {
                    // Try LC correction
                    auto lc_result = __wang_leading_coeff(
                        g_pos, uni_factors, eval, main_v, uni_fac.content);
                    std::cout << "  LC correction: "
                              << (lc_result.success ? "SUCCESS" : "FAILED") << "\n";

                    if (lc_result.success) {
                        std::cout << "  f_scaled = " << lc_result.f_scaled << "\n";
                        std::cout << "  scaled factors:\n";
                        for (size_t i = 0; i < lc_result.scaled_factors.size(); ++i) {
                            std::cout << "    v[" << i << "]: ";
                            for (size_t t = 0; t < lc_result.scaled_factors[i].size(); ++t) {
                                if (t > 0) std::cout << " + ";
                                std::cout << lc_result.scaled_factors[i][t].second
                                          << "*x^" << lc_result.scaled_factors[i][t].first.deg();
                            }
                            std::cout << "\n";
                        }
                        std::cout << "  lc_targets:\n";
                        for (size_t i = 0; i < lc_result.lc_targets.size(); ++i)
                            std::cout << "    tau[" << i << "] = " << lc_result.lc_targets[i] << "\n";

                        // Try MTSHL lift
                        auto mv_factors = __mtshl_lift(
                            lc_result.f_scaled, lc_result.scaled_factors,
                            lc_result.lc_targets, eval, main_v, mtshl_p);

                        std::cout << "  MTSHL lift: "
                                  << (mv_factors.empty() ? "FAILED" : "SUCCESS")
                                  << ", returned " << mv_factors.size() << " factors\n";
                        for (size_t i = 0; i < mv_factors.size(); ++i)
                            std::cout << "    H[" << i << "] = " << mv_factors[i] << "\n";

                        if (!mv_factors.empty()) {
                            // Simulate trial division
                            std::cout << "\n  Trial division (Zassenhaus recombination):\n";

                            auto normalize_factor = [](PolyZZ& h) {
                                h = pp(h);
                                if (!h.empty() && h.front().second < 0)
                                    for (auto& term : h.data())
                                        term.second = -term.second;
                            };

                            std::vector<PolyZZ> normed(mv_factors.size());
                            std::vector<size_t> active_idx;
                            for (size_t fi = 0; fi < mv_factors.size(); ++fi) {
                                normed[fi] = mv_factors[fi];
                                normalize_factor(normed[fi]);
                                std::cout << "    normed[" << fi << "] = " << normed[fi] << "\n";
                                if (!normed[fi].empty() && !is_number(normed[fi]))
                                    active_idx.push_back(fi);
                            }
                            std::cout << "    active_idx count: " << active_idx.size() << "\n";

                            // Try s=1 trial divisions
                            PolyZZ g_remaining = g_pos;
                            std::cout << "\n    s=1 trial divisions:\n";
                            for (size_t ai = 0; ai < active_idx.size(); ++ai) {
                                size_t fi = active_idx[ai];
                                PolyZZ q, rem;
                                pair_vec_div(q.data(), rem.data(),
                                             g_remaining.data(), normed[fi].data(), g_pos.comp());
                                std::cout << "      normed[" << fi << "] divides g_remaining? "
                                          << (rem.empty() && !q.empty() ? "YES" : "NO") << "\n";
                                if (!rem.empty()) {
                                    // Show rem size for diagnosis
                                    std::cout << "        remainder has " << rem.size() << " terms\n";
                                }
                            }

                            // Try s=2 trial divisions
                            if (active_idx.size() >= 4) {
                                std::cout << "\n    s=2 trial divisions (products of pairs):\n";
                                for (size_t i = 0; i < active_idx.size(); ++i) {
                                    for (size_t j = i+1; j < active_idx.size(); ++j) {
                                        PolyZZ prod = normed[active_idx[i]] * normed[active_idx[j]];
                                        prod.normalization();
                                        normalize_factor(prod);
                                        PolyZZ q, rem;
                                        pair_vec_div(q.data(), rem.data(),
                                                     g_remaining.data(), prod.data(), g_pos.comp());
                                        std::cout << "      normed[" << active_idx[i] << "]*normed["
                                                  << active_idx[j] << "] divides? "
                                                  << (rem.empty() && !q.empty() ? "YES" : "NO") << "\n";
                                        if (rem.empty() && !q.empty()) {
                                            std::cout << "        prod = " << prod << "\n";
                                            std::cout << "        quot = " << q << "\n";
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                found_points++;
            }
        }
    }

    // ============================================================
    // Step 3: Call factorize(f) and print full result
    // ============================================================
    std::cout << "\n========================================\n";
    std::cout << "Step 3: factorize(f) - full result\n";
    std::cout << "========================================\n";

    auto fac = factorize(f_lex);
    std::cout << "content = " << fac.content << "\n";
    std::cout << "# factors = " << fac.factors.size() << "\n";
    for (size_t i = 0; i < fac.factors.size(); ++i)
        std::cout << "  factor[" << i << "] (mult=" << fac.factors[i].second << "): "
                  << fac.factors[i].first << "\n";

    // Verify
    PolyZZ check;
    check.data().push_back({{}, fac.content});
    check.normalization();
    for (auto& [fi, ei] : fac.factors)
        for (uint64_t e = 0; e < ei; ++e) {
            check = check * fi;
            check.normalization();
        }
    std::cout << "\nVerification: product == f? " << (check == f_lex ? "YES" : "NO") << "\n";

    // Check if any factor is reducible
    std::cout << "\nChecking reducibility of each factor:\n";
    for (size_t i = 0; i < fac.factors.size(); ++i) {
        auto& fi = fac.factors[i].first;
        auto fi_vars = get_variables(fi);
        if (fi_vars.size() > 1) {
            // Check if fi = f2*f4 or similar
            {
                PolyZZ q, rem;
                pair_vec_div(q.data(), rem.data(),
                             fi.data(), f2_lex.data(), fi.comp());
                if (rem.empty() && !q.empty()) {
                    std::cout << "  factor[" << i << "] is divisible by f2! quotient = " << q << "\n";
                }
            }
            {
                PolyZZ q, rem;
                pair_vec_div(q.data(), rem.data(),
                             fi.data(), f3_lex.data(), fi.comp());
                if (rem.empty() && !q.empty()) {
                    std::cout << "  factor[" << i << "] is divisible by f3! quotient = " << q << "\n";
                }
            }
            {
                PolyZZ q, rem;
                pair_vec_div(q.data(), rem.data(),
                             fi.data(), f4_lex.data(), fi.comp());
                if (rem.empty() && !q.empty()) {
                    std::cout << "  factor[" << i << "] is divisible by f4! quotient = " << q << "\n";
                }
            }

            // Also try factorize on this factor alone
            std::cout << "  factorize(factor[" << i << "]):\n";
            auto sub_fac = factorize(fi);
            std::cout << "    content = " << sub_fac.content << "\n";
            std::cout << "    # sub-factors = " << sub_fac.factors.size() << "\n";
            for (size_t j = 0; j < sub_fac.factors.size(); ++j)
                std::cout << "    sub[" << j << "] (mult=" << sub_fac.factors[j].second << "): "
                          << sub_fac.factors[j].first << "\n";
        }
    }

    return 0;
}
