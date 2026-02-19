/**
 * @file bench_comparative.cc
 * @brief CLPoly vs FLINT (multivariate) and CLPoly vs NTL (univariate) benchmarks.
 *
 * Build with: make bench-comparative
 * Uses release flags (-O3 -DNDEBUG) for meaningful timings.
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mpoly_factor.h>
#include "crosscheck_flint.hh"
#include "crosscheck_ntl.hh"
#include "bench_utils.hh"

using namespace clpoly;

// ---- Memory measurement helpers for FLINT / NTL ----

// Measure per-polynomial memory for fmpz_mpoly (shared context, N copies)
inline double measure_flint_mpoly_kb(int N, const polynomial_ZZ& p,
                                     const std::vector<variable>& vars)
{
    slong nvars = static_cast<slong>(vars.size());
    std::map<variable, slong> var_idx;
    for (slong i = 0; i < nvars; ++i)
        var_idx[vars[i]] = i;

    // Build reference polynomial in shared context
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_ctx_init(ctx, nvars, ORD_DEGLEX);
    fmpz_mpoly_t ref;
    fmpz_mpoly_init(ref, ctx);

    fmpz_t coeff;
    fmpz_init(coeff);
    std::vector<ulong> exp(nvars, 0);
    for (auto& term : p) {
        std::fill(exp.begin(), exp.end(), 0);
        for (auto& ve : term.first) {
            auto it = var_idx.find(ve.first);
            if (it != var_idx.end())
                exp[it->second] = static_cast<ulong>(ve.second);
        }
        crosscheck::clpoly_zz_to_fmpz(coeff, term.second);
        fmpz_mpoly_push_term_fmpz_ui(ref, coeff, exp.data(), ctx);
    }
    fmpz_clear(coeff);
    fmpz_mpoly_sort_terms(ref, ctx);

    // Measure N copies
    long mem0 = get_heap_bytes();
    auto* copies = new fmpz_mpoly_struct[N];
    for (int i = 0; i < N; ++i) {
        fmpz_mpoly_init(&copies[i], ctx);
        fmpz_mpoly_set(&copies[i], ref, ctx);
    }
    long mem1 = get_heap_bytes();
    double kb = (mem1 - mem0) / (1024.0 * N);

    for (int i = 0; i < N; ++i)
        fmpz_mpoly_clear(&copies[i], ctx);
    delete[] copies;
    fmpz_mpoly_clear(ref, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return kb;
}

// Measure per-polynomial memory for fmpz_poly (univariate, N copies)
inline double measure_flint_upoly_kb(int N, const upolynomial_ZZ& p)
{
    auto fp = crosscheck::clpoly_upoly_to_flint(p);

    long mem0 = get_heap_bytes();
    auto* copies = new fmpz_poly_struct[N];
    for (int i = 0; i < N; ++i) {
        fmpz_poly_init(&copies[i]);
        fmpz_poly_set(&copies[i], fp.poly);
    }
    long mem1 = get_heap_bytes();
    double kb = (mem1 - mem0) / (1024.0 * N);

    for (int i = 0; i < N; ++i)
        fmpz_poly_clear(&copies[i]);
    delete[] copies;
    return kb;
}

// Measure per-polynomial memory for NTL::ZZX (N copies)
inline double measure_ntl_zzx_kb(int N, const NTL::ZZX& p)
{
    long mem0 = get_heap_bytes();
    std::vector<NTL::ZZX> copies;
    copies.reserve(N);
    for (int i = 0; i < N; ++i)
        copies.push_back(p);
    long mem1 = get_heap_bytes();
    return (mem1 - mem0) / (1024.0 * N);
}

int main() {
    print_sysinfo();
    variable x("x"), y("y");

    // ================================================================
    // FLINT comparison: Multivariate ZZ
    // ================================================================
    bench_cmp_header("Multivariate ZZ: CLPoly vs FLINT", "FLINT");
    {
        std::vector<variable> vars = {x, y};

        // -- add --
        {
            auto a = random_polynomial<ZZ>(vars, 20, 15, {-50, 50});
            auto b = random_polynomial<ZZ>(vars, 20, 15, {-50, 50});
            auto fa = crosscheck::clpoly_to_flint(a, vars);
            auto fb = crosscheck::clpoly_to_flint(b, vars);

            BENCH_CMP("add  deg20 len15", 100,
                { volatile auto r = a + b; (void)r; },
                {
                    crosscheck::FlintPoly fr(vars);
                    fmpz_mpoly_add(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }
        {
            auto a = random_polynomial<ZZ>(vars, 40, 30, {-50, 50});
            auto b = random_polynomial<ZZ>(vars, 40, 30, {-50, 50});
            auto fa = crosscheck::clpoly_to_flint(a, vars);
            auto fb = crosscheck::clpoly_to_flint(b, vars);

            BENCH_CMP("add  deg40 len30", 100,
                { volatile auto r = a + b; (void)r; },
                {
                    crosscheck::FlintPoly fr(vars);
                    fmpz_mpoly_add(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }

        // -- mul --
        {
            auto a = random_polynomial<ZZ>(vars, 10, 8, {-50, 50});
            auto b = random_polynomial<ZZ>(vars, 10, 8, {-50, 50});
            auto fa = crosscheck::clpoly_to_flint(a, vars);
            auto fb = crosscheck::clpoly_to_flint(b, vars);

            BENCH_CMP("mul  deg10 len8", 50,
                { volatile auto r = a * b; (void)r; },
                {
                    crosscheck::FlintPoly fr(vars);
                    fmpz_mpoly_mul(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }
        {
            auto a = random_polynomial<ZZ>(vars, 20, 15, {-50, 50});
            auto b = random_polynomial<ZZ>(vars, 20, 15, {-50, 50});
            auto fa = crosscheck::clpoly_to_flint(a, vars);
            auto fb = crosscheck::clpoly_to_flint(b, vars);

            BENCH_CMP("mul  deg20 len15", 20,
                { volatile auto r = a * b; (void)r; },
                {
                    crosscheck::FlintPoly fr(vars);
                    fmpz_mpoly_mul(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }

        // -- gcd --
        {
            auto common = random_polynomial<ZZ>(vars, 4, 5, {-10, 10});
            auto a = random_polynomial<ZZ>(vars, 8, 8, {-10, 10}) * common;
            auto b = random_polynomial<ZZ>(vars, 8, 8, {-10, 10}) * common;
            auto va = crosscheck::collect_vars(a, b);
            auto fa = crosscheck::clpoly_to_flint(a, va);
            auto fb = crosscheck::clpoly_to_flint(b, va);

            BENCH_CMP("gcd  deg8+common4", 10,
                { volatile auto r = gcd(a, b); (void)r; },
                {
                    crosscheck::FlintPoly fr(va);
                    fmpz_mpoly_gcd(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }
        {
            auto common = random_polynomial<ZZ>(vars, 8, 8, {-10, 10});
            auto a = random_polynomial<ZZ>(vars, 15, 12, {-10, 10}) * common;
            auto b = random_polynomial<ZZ>(vars, 15, 12, {-10, 10}) * common;
            auto va = crosscheck::collect_vars(a, b);
            auto fa = crosscheck::clpoly_to_flint(a, va);
            auto fb = crosscheck::clpoly_to_flint(b, va);

            BENCH_CMP("gcd  deg15+common8", 5,
                { volatile auto r = gcd(a, b); (void)r; },
                {
                    crosscheck::FlintPoly fr(va);
                    fmpz_mpoly_gcd(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }

        // -- pow --
        {
            auto base = random_polynomial<ZZ>(vars, 3, 4, {-5, 5});
            auto fbase = crosscheck::clpoly_to_flint(base, vars);

            BENCH_CMP("pow  deg3^5", 20,
                { volatile auto r = pow(base, 5); (void)r; },
                {
                    crosscheck::FlintPoly fr(vars);
                    fmpz_mpoly_pow_ui(fr.poly, fbase.poly, 5, fbase.ctx);
                }
            );
            BENCH_CMP("pow  deg3^10", 5,
                { volatile auto r = pow(base, 10); (void)r; },
                {
                    crosscheck::FlintPoly fr(vars);
                    fmpz_mpoly_pow_ui(fr.poly, fbase.poly, 10, fbase.ctx);
                }
            );
        }
    }

    // ================================================================
    // FLINT comparison: QQ mul
    // ================================================================
    bench_cmp_header("Multivariate QQ: CLPoly vs FLINT", "FLINT");
    {
        std::vector<variable> vars = {x, y};
        {
            auto a = random_polynomial_QQ(vars, 10, 8, {-20, 20}, 10);
            auto b = random_polynomial_QQ(vars, 10, 8, {-20, 20}, 10);
            auto fa = crosscheck::clpoly_qq_to_flint(a, vars);
            auto fb = crosscheck::clpoly_qq_to_flint(b, vars);

            BENCH_CMP("mul  deg10", 50,
                { volatile auto r = a * b; (void)r; },
                {
                    crosscheck::FlintQPoly fr(vars);
                    fmpq_mpoly_mul(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }
        {
            auto a = random_polynomial_QQ(vars, 20, 15, {-20, 20}, 10);
            auto b = random_polynomial_QQ(vars, 20, 15, {-20, 20}, 10);
            auto fa = crosscheck::clpoly_qq_to_flint(a, vars);
            auto fb = crosscheck::clpoly_qq_to_flint(b, vars);

            BENCH_CMP("mul  deg20", 10,
                { volatile auto r = a * b; (void)r; },
                {
                    crosscheck::FlintQPoly fr(vars);
                    fmpq_mpoly_mul(fr.poly, fa.poly, fb.poly, fa.ctx);
                }
            );
        }
    }

    // ================================================================
    // NTL comparison: Univariate ZZ
    // ================================================================
    bench_cmp_header("Univariate ZZ: CLPoly vs NTL", "NTL");
    {
        // -- add --
        {
            auto a = random_upolynomial<ZZ>(100, 80, {-100, 100});
            auto b = random_upolynomial<ZZ>(100, 80, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            auto nb = crosscheck::clpoly_upoly_to_ntl(b);

            BENCH_CMP("add  deg100", 100,
                { volatile auto r = a + b; (void)r; },
                { volatile auto r = na + nb; (void)r; }
            );
        }
        {
            auto a = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto b = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            auto nb = crosscheck::clpoly_upoly_to_ntl(b);

            BENCH_CMP("add  deg500", 100,
                { volatile auto r = a + b; (void)r; },
                { volatile auto r = na + nb; (void)r; }
            );
        }

        // -- mul --
        {
            auto a = random_upolynomial<ZZ>(100, 80, {-100, 100});
            auto b = random_upolynomial<ZZ>(100, 80, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            auto nb = crosscheck::clpoly_upoly_to_ntl(b);

            BENCH_CMP("mul  deg100", 50,
                { volatile auto r = a * b; (void)r; },
                { volatile auto r = na * nb; (void)r; }
            );
        }
        {
            auto a = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto b = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            auto nb = crosscheck::clpoly_upoly_to_ntl(b);

            BENCH_CMP("mul  deg500", 10,
                { volatile auto r = a * b; (void)r; },
                { volatile auto r = na * nb; (void)r; }
            );
        }

        // -- gcd --
        {
            auto common = random_upolynomial<ZZ>(25, 20, {-20, 20});
            auto a = random_upolynomial<ZZ>(50, 40, {-20, 20}) * common;
            auto b = random_upolynomial<ZZ>(50, 40, {-20, 20}) * common;
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            auto nb = crosscheck::clpoly_upoly_to_ntl(b);

            BENCH_CMP("gcd  deg50+common25", 10,
                { volatile auto r = polynomial_GCD(a, b); (void)r; },
                { volatile auto r = NTL::GCD(na, nb); (void)r; }
            );
        }
        {
            auto common = random_upolynomial<ZZ>(100, 80, {-20, 20});
            auto a = random_upolynomial<ZZ>(200, 150, {-20, 20}) * common;
            auto b = random_upolynomial<ZZ>(200, 150, {-20, 20}) * common;
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            auto nb = crosscheck::clpoly_upoly_to_ntl(b);

            BENCH_CMP("gcd  deg200+common100", 3,
                { volatile auto r = polynomial_GCD(a, b); (void)r; },
                { volatile auto r = NTL::GCD(na, nb); (void)r; }
            );
        }

        // -- derivative --
        {
            auto a = random_upolynomial<ZZ>(100, 80, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);

            BENCH_CMP("deriv  deg100", 100,
                { volatile auto r = derivative(a); (void)r; },
                { volatile auto r = crosscheck::ntl_diff(na); (void)r; }
            );
        }
        {
            auto a = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);

            BENCH_CMP("deriv  deg500", 100,
                { volatile auto r = derivative(a); (void)r; },
                { volatile auto r = crosscheck::ntl_diff(na); (void)r; }
            );
        }

        // -- eval --
        {
            auto a = random_upolynomial<ZZ>(100, 80, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            ZZ pt(42);
            NTL::ZZ npt = crosscheck::clpoly_zz_to_ntl(pt);

            BENCH_CMP("eval  deg100", 100,
                { volatile auto r = assign(a, pt); (void)r; },
                { volatile auto r = crosscheck::ntl_eval(na, npt); (void)r; }
            );
        }
        {
            auto a = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto na = crosscheck::clpoly_upoly_to_ntl(a);
            ZZ pt(42);
            NTL::ZZ npt = crosscheck::clpoly_zz_to_ntl(pt);

            BENCH_CMP("eval  deg500", 50,
                { volatile auto r = assign(a, pt); (void)r; },
                { volatile auto r = crosscheck::ntl_eval(na, npt); (void)r; }
            );
        }
    }

    // ================================================================
    // Factorization: CLPoly vs FLINT vs NTL (three-way)
    // ================================================================
    bench_cmp_header("Factorization: CLPoly vs FLINT", "FLINT");
    {
        auto f1 = random_upolynomial<ZZ>(5, 4, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(5, 4, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(5, 4, {-10, 10});
        auto f4 = random_upolynomial<ZZ>(6, 5, {-10, 10});
        auto f5 = random_upolynomial<ZZ>(8, 6, {-10, 10});

        auto p_s = f1 * f2 * f3;
        auto p_m = f1 * f2 * f3 * f4;
        auto p_l = f1 * f2 * f3 * f4 * f5;

        BENCH_CMP("factor  ~deg15 (3 fac)", 5,
            { volatile auto r = factorize(p_s); (void)r; },
            { volatile auto r = crosscheck::flint_factor_upoly(p_s); (void)r; }
        );
        BENCH_CMP("factor  ~deg21 (4 fac)", 3,
            { volatile auto r = factorize(p_m); (void)r; },
            { volatile auto r = crosscheck::flint_factor_upoly(p_m); (void)r; }
        );
        BENCH_CMP("factor  ~deg29 (5 fac)", 3,
            { volatile auto r = factorize(p_l); (void)r; },
            { volatile auto r = crosscheck::flint_factor_upoly(p_l); (void)r; }
        );
    }

    bench_cmp_header("Factorization: CLPoly vs NTL", "NTL");
    {
        auto f1 = random_upolynomial<ZZ>(5, 4, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(5, 4, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(5, 4, {-10, 10});
        auto f4 = random_upolynomial<ZZ>(6, 5, {-10, 10});
        auto f5 = random_upolynomial<ZZ>(8, 6, {-10, 10});

        auto p_s = f1 * f2 * f3;
        auto p_m = f1 * f2 * f3 * f4;
        auto p_l = f1 * f2 * f3 * f4 * f5;
        auto ns = crosscheck::clpoly_upoly_to_ntl(p_s);
        auto nm = crosscheck::clpoly_upoly_to_ntl(p_m);
        auto nl = crosscheck::clpoly_upoly_to_ntl(p_l);

        BENCH_CMP("factor  ~deg15 (3 fac)", 5,
            { volatile auto r = factorize(p_s); (void)r; },
            { volatile auto r = crosscheck::ntl_factor(ns); (void)r; }
        );
        BENCH_CMP("factor  ~deg21 (4 fac)", 3,
            { volatile auto r = factorize(p_m); (void)r; },
            { volatile auto r = crosscheck::ntl_factor(nm); (void)r; }
        );
        BENCH_CMP("factor  ~deg29 (5 fac)", 3,
            { volatile auto r = factorize(p_l); (void)r; },
            { volatile auto r = crosscheck::ntl_factor(nl); (void)r; }
        );
    }

    // ================================================================
    // Multivariate Factorization: CLPoly vs FLINT
    // ================================================================
    bench_cmp_header("Multivariate Factorize: CLPoly vs FLINT", "FLINT");
    {
        using PolyLex = polynomial_<ZZ, lex>;
        auto to_lex = [](const polynomial_ZZ& p) {
            PolyLex r; poly_convert(p, r); return r;
        };

        std::vector<variable> vars = {x, y};

        // Known: x^2 - y^2
        {
            auto f_lex = to_lex(pow(x,2) - pow(y,2));
            auto ff = crosscheck::clpoly_to_flint(pow(x,2) - pow(y,2), vars);

            BENCH_CMP("factor  x^2-y^2", 10,
                { volatile auto r = factorize(f_lex); (void)r; },
                {
                    fmpz_mpoly_factor_t fac;
                    fmpz_mpoly_factor_init(fac, ff.ctx);
                    fmpz_mpoly_factor(fac, ff.poly, ff.ctx);
                    fmpz_mpoly_factor_clear(fac, ff.ctx);
                }
            );
        }

        // Bivariate: product of two random deg3 factors
        {
            auto f1_gr = random_polynomial<ZZ>(vars, 3, 4, {-5, 5});
            auto f2_gr = random_polynomial<ZZ>(vars, 3, 4, {-5, 5});
            auto f_gr = f1_gr * f2_gr;
            auto f_lex = to_lex(f_gr);
            auto ff = crosscheck::clpoly_to_flint(f_gr, vars);

            BENCH_CMP("factor  bivar deg3*deg3", 5,
                { volatile auto r = factorize(f_lex); (void)r; },
                {
                    fmpz_mpoly_factor_t fac;
                    fmpz_mpoly_factor_init(fac, ff.ctx);
                    fmpz_mpoly_factor(fac, ff.poly, ff.ctx);
                    fmpz_mpoly_factor_clear(fac, ff.ctx);
                }
            );
        }

        // Bivariate: product of two random deg5 factors
        {
            auto f1_gr = random_polynomial<ZZ>(vars, 5, 6, {-5, 5});
            auto f2_gr = random_polynomial<ZZ>(vars, 5, 6, {-5, 5});
            auto f_gr = f1_gr * f2_gr;
            auto f_lex = to_lex(f_gr);
            auto ff = crosscheck::clpoly_to_flint(f_gr, vars);

            BENCH_CMP("factor  bivar deg5*deg5", 3,
                { volatile auto r = factorize(f_lex); (void)r; },
                {
                    fmpz_mpoly_factor_t fac;
                    fmpz_mpoly_factor_init(fac, ff.ctx);
                    fmpz_mpoly_factor(fac, ff.poly, ff.ctx);
                    fmpz_mpoly_factor_clear(fac, ff.ctx);
                }
            );
        }

        // Trivariate: known (x+y+z)(x-y+z)
        {
            variable z("z");
            std::vector<variable> vars3 = {x, y, z};
            auto f_gr = pow(x,2) + 2*x*z + pow(z,2) - pow(y,2);
            auto f_lex = to_lex(f_gr);
            auto ff = crosscheck::clpoly_to_flint(f_gr, vars3);

            BENCH_CMP("factor  trivar known", 5,
                { volatile auto r = factorize(f_lex); (void)r; },
                {
                    fmpz_mpoly_factor_t fac;
                    fmpz_mpoly_factor_init(fac, ff.ctx);
                    fmpz_mpoly_factor(fac, ff.poly, ff.ctx);
                    fmpz_mpoly_factor_clear(fac, ff.ctx);
                }
            );
        }
    }

    // ================================================================
    // Memory Comparison: per-polynomial representation size
    // ================================================================
    const int MEM_N = 500;

    // -- Multivariate ZZ: CLPoly vs FLINT --
    bench_header("Memory: CLPoly polynomial_ZZ vs FLINT fmpz_mpoly");
    mem_cmp_header("FLINT");
    {
        std::vector<variable> vars = {x, y};
        {
            auto p = random_polynomial<ZZ>(vars, 20, 15, {-50, 50});
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_ZZ(p); });
            double fl_kb = measure_flint_mpoly_kb(MEM_N, p, vars);
            mem_cmp_row("mpoly_ZZ deg20 len15", p.size(), cl_kb, fl_kb);
        }
        {
            auto p = random_polynomial<ZZ>(vars, 40, 30, {-50, 50});
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_ZZ(p); });
            double fl_kb = measure_flint_mpoly_kb(MEM_N, p, vars);
            mem_cmp_row("mpoly_ZZ deg40 len30", p.size(), cl_kb, fl_kb);
        }
        {
            auto p = random_polynomial<ZZ>(vars, 80, 60, {-50, 50});
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_ZZ(p); });
            double fl_kb = measure_flint_mpoly_kb(MEM_N, p, vars);
            mem_cmp_row("mpoly_ZZ deg80 len60", p.size(), cl_kb, fl_kb);
        }
    }

    // -- Univariate ZZ: CLPoly vs FLINT --
    bench_header("Memory: CLPoly upolynomial_ZZ vs FLINT fmpz_poly");
    mem_cmp_header("FLINT");
    {
        {
            auto p = random_upolynomial<ZZ>(50, 40, {-100, 100});
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
            double fl_kb = measure_flint_upoly_kb(MEM_N, p);
            mem_cmp_row("upoly_ZZ deg50", p.size(), cl_kb, fl_kb);
        }
        {
            auto p = random_upolynomial<ZZ>(200, 150, {-100, 100});
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
            double fl_kb = measure_flint_upoly_kb(MEM_N, p);
            mem_cmp_row("upoly_ZZ deg200", p.size(), cl_kb, fl_kb);
        }
        {
            auto p = random_upolynomial<ZZ>(500, 400, {-100, 100});
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
            double fl_kb = measure_flint_upoly_kb(MEM_N, p);
            mem_cmp_row("upoly_ZZ deg500", p.size(), cl_kb, fl_kb);
        }
    }

    // -- Univariate ZZ: CLPoly vs NTL --
    bench_header("Memory: CLPoly upolynomial_ZZ vs NTL ZZX");
    mem_cmp_header("NTL");
    {
        {
            auto p = random_upolynomial<ZZ>(50, 40, {-100, 100});
            auto np = crosscheck::clpoly_upoly_to_ntl(p);
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
            double ntl_kb = measure_ntl_zzx_kb(MEM_N, np);
            mem_cmp_row("upoly_ZZ deg50", p.size(), cl_kb, ntl_kb);
        }
        {
            auto p = random_upolynomial<ZZ>(200, 150, {-100, 100});
            auto np = crosscheck::clpoly_upoly_to_ntl(p);
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
            double ntl_kb = measure_ntl_zzx_kb(MEM_N, np);
            mem_cmp_row("upoly_ZZ deg200", p.size(), cl_kb, ntl_kb);
        }
        {
            auto p = random_upolynomial<ZZ>(500, 400, {-100, 100});
            auto np = crosscheck::clpoly_upoly_to_ntl(p);
            double cl_kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
            double ntl_kb = measure_ntl_zzx_kb(MEM_N, np);
            mem_cmp_row("upoly_ZZ deg500", p.size(), cl_kb, ntl_kb);
        }
    }

    print_process_memory();
    std::cout << "\nDone. Total benchmarks: " << _bench_results.size() << std::endl;
    return 0;
}
