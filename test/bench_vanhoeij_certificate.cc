// Fixed-prime benchmark: every integer quadratic is irreducible (x^2+c,
// c>0), but splits into two distinct linear factors modulo 97.  Pass "hard"
// for two Eisenstein irreducibles with r/2 modular factors each, or
// "nonmonic" for primitive positive nonmonic quadratics with the same roots.
#include <clpoly/polynomial_factorize_univar.hh>
#include <chrono>
#include <iostream>
#include <string>
#ifdef CLPOLY_FLINT_REFERENCE
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include "crosscheck_flint.hh"
#endif

using namespace clpoly;

int main(int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : "quadratic";
    if (mode != "quadratic" && mode != "hard" && mode != "nonmonic") return 8;
    const std::string selected_r = argc > 2 ? argv[2] : "all";
    if (selected_r != "all" && selected_r != "24" && selected_r != "32") return 9;
    const bool hard = mode == "hard";
    const bool nonmonic = mode == "nonmonic";
    for (int r : {24, 32})
    {
        if (selected_r != "all" && selected_r != std::to_string(r)) continue;
        upolynomial_<ZZ> f;
        f.push_back({umonomial(0), ZZ(1)});
        std::vector<upolynomial_<ZZ>> expected;
        if (hard)
        {
            // Each product has even roots.  Adding 2*97 makes its constant
            // coefficient 2 mod 4, so Eisenstein at 2 proves irreducibility.
            // Modulo 97 the roots remain distinct: two true factors, r
            // modular factors, requiring large subsets in Zassenhaus.
            for (int group = 0; group < 2; ++group)
            {
                auto g = f;
                for (int i = 1; i <= r / 2; ++i)
                {
                    upolynomial_<ZZ> linear;
                    linear.push_back({umonomial(1), ZZ(1)});
                    linear.push_back({umonomial(0), ZZ(-2*(i + group*r/2))});
                    g = g * linear;
                    g.normalization();
                }
                upolynomial_<ZZ> offset;
                offset.push_back({umonomial(0), ZZ(194)});
                g = g + offset;
                g.normalization();
                expected.push_back(g);
            }
            for (const auto& g : expected)
            {
                f = f * g;
                f.normalization();
            }
        }
        else for (int root = 1; root <= r / 2; ++root)
        {
            upolynomial_<ZZ> g;
            const int leading = nonmonic ? 2 * root + 1 : 1;
            g.push_back({umonomial(2), ZZ(leading)});
            g.push_back({umonomial(0), ZZ(97 - leading * root * root % 97)});
            if (nonmonic) g = __upoly_primitive(std::move(g)).second;
            expected.push_back(g);
            f = f * g;
            f.normalization();
        }
        auto modular = __factor_Zp(polynomial_mod(f, 97));
        std::vector<upolynomial_<Zp>> factors;
        for (const auto& entry : modular.second) factors.push_back(entry.first);
        if ((int)factors.size() != r) return 1;
#ifdef CLPOLY_FLINT_REFERENCE
        // Feed the same prime and the same modular factors to the reference
        // van Hoeij entry, bypassing its independent prime-selection policy.
        auto reference_input = crosscheck::clpoly_upoly_to_flint(f);
        nmod_poly_factor_t reference_factors;
        nmod_poly_factor_init(reference_factors);
        nmod_poly_t local;
        nmod_poly_init(local, 97);
        for (const auto& g : factors)
        {
            nmod_poly_zero(local);
            for (const auto& term : g)
                nmod_poly_set_coeff_ui(local, term.first.deg(), term.second.number());
            nmod_poly_factor_insert(reference_factors, local, 1);
        }
        nmod_poly_clear(local);
#endif
        for (int repeat = 0; repeat < 3; ++repeat)
        {
#ifdef CLPOLY_PROFILE
            reset_factorize_profile();
#endif
#ifdef CLPOLY_PHASE_BENCH
            __bench_phases = {};
#endif
            auto start = std::chrono::steady_clock::now();
            auto result = __lll_factorize(f, factors, 97);
            auto end = std::chrono::steady_clock::now();
            if (result.size() != expected.size()) return 2;
            for (const auto& g : expected)
                if (std::count(result.begin(), result.end(), g) != 1) return 3;
            std::cout << "r=" << r << " repeat=" << repeat << " ms="
                      << std::chrono::duration<double, std::milli>(end-start).count()
                      << '\n' << std::flush;
#ifdef CLPOLY_PHASE_BENCH
            std::cout << "phases r=" << r << " repeat=" << repeat
                      << " hensel_ns=" << __bench_phases.hensel_ns
                      << " lll_ns=" << __bench_phases.lll_ns
                      << " inner_zass_ns=" << __bench_phases.inner_zass_ns
                      << " outer_zass_ns=" << __bench_phases.outer_zass_ns
                      << " hensel_calls=" << __bench_phases.hensel_calls
                      << " lll_calls=" << __bench_phases.lll_calls
                      << " inner_zass_calls=" << __bench_phases.inner_zass_calls
                      << " outer_zass_calls=" << __bench_phases.outer_zass_calls
                      << '\n';
#endif
#ifdef CLPOLY_CERTIFICATE_BENCH
            if (__g_profile.certified_recombinations != 1 ||
                __g_profile.zassenhaus_fallbacks != 0) return 4;
            std::cout << "certified=" << __g_profile.certified_recombinations
                      << " fallback=" << __g_profile.zassenhaus_fallbacks
                      << " lifts=" << 1 + __g_profile.phase2_triggers
                      << " lll=" << __g_profile.lattice_rounds
                      << " columns=" << __g_profile.lattice_columns
                      << '\n';
#endif
#ifdef CLPOLY_FLINT_REFERENCE
            fmpz_poly_factor_t reference_result;
            fmpz_poly_factor_init(reference_result);
            auto reference_start = std::chrono::steady_clock::now();
            fmpz_poly_factor_van_hoeij(reference_result, reference_factors,
                                     reference_input.poly, 1, 97);
            auto reference_end = std::chrono::steady_clock::now();
            if (reference_result->num != static_cast<slong>(expected.size())) return 5;
            std::vector<upolynomial_<ZZ>> reference_output;
            for (slong i = 0; i < reference_result->num; ++i)
            {
                if (reference_result->exp[i] != 1) return 6;
                reference_output.push_back(crosscheck::flint_to_clpoly_upoly(reference_result->p + i));
            }
            for (const auto& g : expected)
                if (std::count(reference_output.begin(), reference_output.end(), g) != 1) return 7;
            std::cout << "flint=" << FLINT_VERSION << " r=" << r << " repeat=" << repeat
                      << " ms=" << std::chrono::duration<double, std::milli>(reference_end-reference_start).count()
                      << " factors=" << reference_result->num << '\n' << std::flush;
            fmpz_poly_factor_clear(reference_result);
#endif
        }
#ifdef CLPOLY_FLINT_REFERENCE
        nmod_poly_factor_clear(reference_factors);
#endif
    }
}
