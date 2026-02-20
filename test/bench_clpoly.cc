/**
 * @file bench_clpoly.cc
 * @brief Pure CLPoly performance benchmarks covering all major operations.
 *
 * Build with: make bench-clpoly
 * Uses release flags (-O3 -DNDEBUG) for meaningful timings.
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include "bench_utils.hh"

using namespace clpoly;

int main() {
    print_sysinfo();
    variable x("x"), y("y");

    // ================================================================
    // Multivariate polynomial_ZZ (2 variables)
    // ================================================================
    bench_header("Multivariate polynomial_ZZ (2 variables)");

    // Pre-generate test data outside BENCH to measure only the operation
    // Small
    auto mz_a_s = random_polynomial<ZZ>({x, y}, 20, 15, {-50, 50});
    auto mz_b_s = random_polynomial<ZZ>({x, y}, 20, 15, {-50, 50});
    // Medium
    auto mz_a_m = random_polynomial<ZZ>({x, y}, 40, 30, {-50, 50});
    auto mz_b_m = random_polynomial<ZZ>({x, y}, 40, 30, {-50, 50});
    // Large
    auto mz_a_l = random_polynomial<ZZ>({x, y}, 80, 60, {-50, 50});
    auto mz_b_l = random_polynomial<ZZ>({x, y}, 80, 60, {-50, 50});

    // -- add --
    BENCH("add  deg20 len15", 100, {
        volatile auto r = mz_a_s + mz_b_s; (void)r;
    });
    BENCH("add  deg40 len30", 100, {
        volatile auto r = mz_a_m + mz_b_m; (void)r;
    });
    BENCH("add  deg80 len60", 100, {
        volatile auto r = mz_a_l + mz_b_l; (void)r;
    });

    // -- mul --
    auto mz_mul_a_s = random_polynomial<ZZ>({x, y}, 10, 8, {-50, 50});
    auto mz_mul_b_s = random_polynomial<ZZ>({x, y}, 10, 8, {-50, 50});
    auto mz_mul_a_m = random_polynomial<ZZ>({x, y}, 20, 15, {-50, 50});
    auto mz_mul_b_m = random_polynomial<ZZ>({x, y}, 20, 15, {-50, 50});
    auto mz_mul_a_l = random_polynomial<ZZ>({x, y}, 40, 30, {-50, 50});
    auto mz_mul_b_l = random_polynomial<ZZ>({x, y}, 40, 30, {-50, 50});

    BENCH("mul  deg10 len8", 50, {
        volatile auto r = mz_mul_a_s * mz_mul_b_s; (void)r;
    });
    BENCH("mul  deg20 len15", 20, {
        volatile auto r = mz_mul_a_m * mz_mul_b_m; (void)r;
    });
    BENCH("mul  deg40 len30", 5, {
        volatile auto r = mz_mul_a_l * mz_mul_b_l; (void)r;
    });

    // -- pow --
    auto mz_pow_base = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
    BENCH("pow  deg3^5", 20, {
        volatile auto r = pow(mz_pow_base, 5); (void)r;
    });
    BENCH("pow  deg3^10", 5, {
        volatile auto r = pow(mz_pow_base, 10); (void)r;
    });
    BENCH("pow  deg3^15", 3, {
        volatile auto r = pow(mz_pow_base, 15); (void)r;
    });

    // -- gcd --
    auto mz_common_s = random_polynomial<ZZ>({x, y}, 4, 5, {-10, 10});
    auto mz_gcd_f1_s = random_polynomial<ZZ>({x, y}, 8, 8, {-10, 10}) * mz_common_s;
    auto mz_gcd_f2_s = random_polynomial<ZZ>({x, y}, 8, 8, {-10, 10}) * mz_common_s;

    auto mz_common_m = random_polynomial<ZZ>({x, y}, 8, 8, {-10, 10});
    auto mz_gcd_f1_m = random_polynomial<ZZ>({x, y}, 15, 12, {-10, 10}) * mz_common_m;
    auto mz_gcd_f2_m = random_polynomial<ZZ>({x, y}, 15, 12, {-10, 10}) * mz_common_m;

    auto mz_common_l = random_polynomial<ZZ>({x, y}, 12, 12, {-10, 10});
    auto mz_gcd_f1_l = random_polynomial<ZZ>({x, y}, 25, 20, {-10, 10}) * mz_common_l;
    auto mz_gcd_f2_l = random_polynomial<ZZ>({x, y}, 25, 20, {-10, 10}) * mz_common_l;

    BENCH("gcd  deg8+common4", 10, {
        volatile auto r = gcd(mz_gcd_f1_s, mz_gcd_f2_s); (void)r;
    });
    BENCH("gcd  deg15+common8", 5, {
        volatile auto r = gcd(mz_gcd_f1_m, mz_gcd_f2_m); (void)r;
    });
    BENCH("gcd  deg25+common12", 3, {
        volatile auto r = gcd(mz_gcd_f1_l, mz_gcd_f2_l); (void)r;
    });

    // -- derivative --
    BENCH("deriv  deg20", 100, {
        volatile auto r = derivative(mz_a_s, x); (void)r;
    });
    BENCH("deriv  deg40", 100, {
        volatile auto r = derivative(mz_a_m, x); (void)r;
    });
    BENCH("deriv  deg80", 100, {
        volatile auto r = derivative(mz_a_l, x); (void)r;
    });

    // -- prem --
    auto mz_prem_f_s = random_polynomial<ZZ>({x, y}, 8, 8, {-20, 20});
    auto mz_prem_g_s = random_polynomial<ZZ>({x, y}, 5, 5, {-20, 20});
    auto mz_prem_f_m = random_polynomial<ZZ>({x, y}, 15, 12, {-20, 20});
    auto mz_prem_g_m = random_polynomial<ZZ>({x, y}, 8, 8, {-20, 20});
    auto mz_prem_f_l = random_polynomial<ZZ>({x, y}, 25, 20, {-20, 20});
    auto mz_prem_g_l = random_polynomial<ZZ>({x, y}, 12, 12, {-20, 20});

    BENCH("prem  deg8/5", 20, {
        volatile auto r = prem(mz_prem_f_s, mz_prem_g_s, x); (void)r;
    });
    BENCH("prem  deg15/8", 10, {
        volatile auto r = prem(mz_prem_f_m, mz_prem_g_m, x); (void)r;
    });
    BENCH("prem  deg25/12", 5, {
        volatile auto r = prem(mz_prem_f_l, mz_prem_g_l, x); (void)r;
    });

    // -- resultant --
    auto mz_res_f_s = random_polynomial<ZZ>({x, y}, 4, 5, {-10, 10});
    auto mz_res_g_s = random_polynomial<ZZ>({x, y}, 4, 5, {-10, 10});
    auto mz_res_f_m = random_polynomial<ZZ>({x, y}, 6, 8, {-10, 10});
    auto mz_res_g_m = random_polynomial<ZZ>({x, y}, 6, 8, {-10, 10});
    auto mz_res_f_l = random_polynomial<ZZ>({x, y}, 8, 10, {-10, 10});
    auto mz_res_g_l = random_polynomial<ZZ>({x, y}, 8, 10, {-10, 10});

    BENCH("resultant  deg4", 10, {
        volatile auto r = resultant(mz_res_f_s, mz_res_g_s, x); (void)r;
    });
    BENCH("resultant  deg6", 5, {
        volatile auto r = resultant(mz_res_f_m, mz_res_g_m, x); (void)r;
    });
    BENCH("resultant  deg8", 3, {
        volatile auto r = resultant(mz_res_f_l, mz_res_g_l, x); (void)r;
    });

    // ================================================================
    // Multivariate factorize (lex ordering)
    // ================================================================
    bench_header("Multivariate factorize (lex)");
    {
        using PolyLex = polynomial_<ZZ, lex>;
        auto to_lex = [](const polynomial_ZZ& p) {
            PolyLex r; poly_convert(p, r); return r;
        };

        // Bivariate: (f1)(f2), small factors
        auto mfac_f1_s = to_lex(random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5}));
        auto mfac_f2_s = to_lex(random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5}));
        auto mfac_s = mfac_f1_s * mfac_f2_s;
        mfac_s.normalization();

        auto mfac_f1_m = to_lex(random_polynomial<ZZ>({x, y}, 5, 6, {-5, 5}));
        auto mfac_f2_m = to_lex(random_polynomial<ZZ>({x, y}, 5, 6, {-5, 5}));
        auto mfac_m = mfac_f1_m * mfac_f2_m;
        mfac_m.normalization();

        // Known factorable: x^2 - y^2 = (x+y)(x-y)
        auto mfac_known = to_lex(pow(x,2) - pow(y,2));

        BENCH("factorize  x^2-y^2", 10, {
            volatile auto r = factorize(mfac_known); (void)r;
        });
        BENCH("factorize  bivar deg3*deg3", 5, {
            volatile auto r = factorize(mfac_s); (void)r;
        });
        BENCH("factorize  bivar deg5*deg5", 3, {
            volatile auto r = factorize(mfac_m); (void)r;
        });

        // Trivariate
        variable z("z");
        auto mfac_tri_f1 = to_lex(random_polynomial<ZZ>({x, y, z}, 2, 4, {-3, 3}));
        auto mfac_tri_f2 = to_lex(random_polynomial<ZZ>({x, y, z}, 2, 4, {-3, 3}));
        auto mfac_tri = mfac_tri_f1 * mfac_tri_f2;
        mfac_tri.normalization();

        auto mfac_tri_known = to_lex(pow(x,2) + 2*x*z + pow(z,2) - pow(y,2));  // (x+y+z)(x-y+z)

        BENCH("factorize  trivar known", 5, {
            volatile auto r = factorize(mfac_tri_known); (void)r;
        });
        BENCH("factorize  trivar deg2*deg2", 3, {
            volatile auto r = factorize(mfac_tri); (void)r;
        });

        // With multiplicity: (x+y)^2 * (x-y)
        auto mfac_mult = to_lex(pow(x,3) + pow(x,2)*y - x*pow(y,2) - pow(y,3));
        BENCH("factorize  (x+y)^2*(x-y)", 10, {
            volatile auto r = factorize(mfac_mult); (void)r;
        });

        // Classic: x^4 - y^4 = (x-y)(x+y)(x^2+y^2)
        auto mfac_x4y4 = to_lex(pow(x,4) - pow(y,4));
        BENCH("factorize  x^4-y^4", 10, {
            volatile auto r = factorize(mfac_x4y4); (void)r;
        });

        // Classic: SymPy f_1 = (x+yz+20)(xy+z+10)(xz+y+30)
        {
            variable z("z");
            auto sympy_f1 = to_lex(
                (x + y*z + polynomial_ZZ(ZZ(20)))
              * (x*y + z + polynomial_ZZ(ZZ(10)))
              * (x*z + y + polynomial_ZZ(ZZ(30))));
            BENCH("factorize  SymPy f_1 (3var 3fac)", 3, {
                volatile auto r = factorize(sympy_f1); (void)r;
            });
        }

        // Classic: (x+y+z)^3 - (x+y+z) = (x+y+z)(x+y+z-1)(x+y+z+1)
        {
            variable z("z");
            auto s = x + y + z;
            auto cube_minus = to_lex(pow(s, 3) - s);
            BENCH("factorize  (x+y+z)^3-(x+y+z)", 5, {
                volatile auto r = factorize(cube_minus); (void)r;
            });
        }
    }

    // ================================================================
    // Univariate upolynomial_ZZ
    // ================================================================
    bench_header("Univariate upolynomial_ZZ");

    auto uz_a_s = random_upolynomial<ZZ>(50, 40, {-100, 100});
    auto uz_b_s = random_upolynomial<ZZ>(50, 40, {-100, 100});
    auto uz_a_m = random_upolynomial<ZZ>(200, 150, {-100, 100});
    auto uz_b_m = random_upolynomial<ZZ>(200, 150, {-100, 100});
    auto uz_a_l = random_upolynomial<ZZ>(500, 400, {-100, 100});
    auto uz_b_l = random_upolynomial<ZZ>(500, 400, {-100, 100});

    // -- add --
    BENCH("add  deg50", 100, {
        volatile auto r = uz_a_s + uz_b_s; (void)r;
    });
    BENCH("add  deg200", 100, {
        volatile auto r = uz_a_m + uz_b_m; (void)r;
    });
    BENCH("add  deg500", 100, {
        volatile auto r = uz_a_l + uz_b_l; (void)r;
    });

    // -- mul --
    BENCH("mul  deg50", 100, {
        volatile auto r = uz_a_s * uz_b_s; (void)r;
    });
    BENCH("mul  deg200", 20, {
        volatile auto r = uz_a_m * uz_b_m; (void)r;
    });
    BENCH("mul  deg500", 5, {
        volatile auto r = uz_a_l * uz_b_l; (void)r;
    });

    // -- gcd --
    auto uz_common_s = random_upolynomial<ZZ>(15, 12, {-20, 20});
    auto uz_gcd_f1_s = random_upolynomial<ZZ>(30, 25, {-20, 20}) * uz_common_s;
    auto uz_gcd_f2_s = random_upolynomial<ZZ>(30, 25, {-20, 20}) * uz_common_s;

    auto uz_common_m = random_upolynomial<ZZ>(40, 30, {-20, 20});
    auto uz_gcd_f1_m = random_upolynomial<ZZ>(80, 60, {-20, 20}) * uz_common_m;
    auto uz_gcd_f2_m = random_upolynomial<ZZ>(80, 60, {-20, 20}) * uz_common_m;

    auto uz_common_l = random_upolynomial<ZZ>(100, 80, {-20, 20});
    auto uz_gcd_f1_l = random_upolynomial<ZZ>(200, 150, {-20, 20}) * uz_common_l;
    auto uz_gcd_f2_l = random_upolynomial<ZZ>(200, 150, {-20, 20}) * uz_common_l;

    BENCH("gcd  deg30+common15", 10, {
        volatile auto r = polynomial_GCD(uz_gcd_f1_s, uz_gcd_f2_s); (void)r;
    });
    BENCH("gcd  deg80+common40", 5, {
        volatile auto r = polynomial_GCD(uz_gcd_f1_m, uz_gcd_f2_m); (void)r;
    });
    BENCH("gcd  deg200+common100", 3, {
        volatile auto r = polynomial_GCD(uz_gcd_f1_l, uz_gcd_f2_l); (void)r;
    });

    // -- derivative --
    BENCH("deriv  deg50", 100, {
        volatile auto r = derivative(uz_a_s); (void)r;
    });
    BENCH("deriv  deg200", 100, {
        volatile auto r = derivative(uz_a_m); (void)r;
    });
    BENCH("deriv  deg500", 100, {
        volatile auto r = derivative(uz_a_l); (void)r;
    });

    // -- eval (assign) --
    ZZ eval_pt(42);
    BENCH("eval  deg50", 100, {
        volatile auto r = assign(uz_a_s, eval_pt); (void)r;
    });
    BENCH("eval  deg200", 50, {
        volatile auto r = assign(uz_a_m, eval_pt); (void)r;
    });
    BENCH("eval  deg500", 20, {
        volatile auto r = assign(uz_a_l, eval_pt); (void)r;
    });

    // -- content --
    auto uz_cont_s = random_upolynomial<ZZ>(50, 40, {-100, 100});
    auto uz_cont_m = random_upolynomial<ZZ>(200, 150, {-100, 100});
    auto uz_cont_l = random_upolynomial<ZZ>(500, 400, {-100, 100});

    BENCH("cont  deg50", 100, {
        volatile auto r = cont(uz_cont_s); (void)r;
    });
    BENCH("cont  deg200", 50, {
        volatile auto r = cont(uz_cont_m); (void)r;
    });
    BENCH("cont  deg500", 20, {
        volatile auto r = cont(uz_cont_l); (void)r;
    });

    // -- factorize --
    // Build products of small irreducible-ish factors
    auto uf1 = random_upolynomial<ZZ>(5, 4, {-10, 10});
    auto uf2 = random_upolynomial<ZZ>(5, 4, {-10, 10});
    auto uf3 = random_upolynomial<ZZ>(5, 4, {-10, 10});
    auto uf4 = random_upolynomial<ZZ>(6, 5, {-10, 10});
    auto uf5 = random_upolynomial<ZZ>(8, 6, {-10, 10});
    auto uz_fac_s = uf1 * uf2 * uf3;                    // ~deg15, 3 factors
    auto uz_fac_m = uf1 * uf2 * uf3 * uf4;              // ~deg21, 4 factors
    auto uz_fac_l = uf1 * uf2 * uf3 * uf4 * uf5;        // ~deg29, 5 factors

    BENCH("factorize  ~deg15 (3 factors)", 5, {
        volatile auto r = factorize(uz_fac_s); (void)r;
    });
    BENCH("factorize  ~deg21 (4 factors)", 3, {
        volatile auto r = factorize(uz_fac_m); (void)r;
    });
    BENCH("factorize  ~deg29 (5 factors)", 3, {
        volatile auto r = factorize(uz_fac_l); (void)r;
    });

    // -- classic univariate factorization --

    // Wilkinson W(10): (x-1)(x-2)...(x-10)
    {
        upolynomial_ZZ wilk10({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 10; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            wilk10 = wilk10 * lin;
        }
        BENCH("factorize  Wilkinson W(10)", 5, {
            volatile auto r = factorize(wilk10); (void)r;
        });
    }

    // Wilkinson W(15): (x-1)(x-2)...(x-15)
    {
        upolynomial_ZZ wilk15({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 15; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            wilk15 = wilk15 * lin;
        }
        BENCH("factorize  Wilkinson W(15)", 3, {
            volatile auto r = factorize(wilk15); (void)r;
        });
    }

    // Cyclotomic x^15-1 (4 irreducible factors)
    {
        upolynomial_ZZ cyc15({{15, ZZ(1)}, {0, ZZ(-1)}});
        BENCH("factorize  x^15-1 (cyclotomic)", 5, {
            volatile auto r = factorize(cyc15); (void)r;
        });
    }

    // Cyclotomic x^24-1 (8 irreducible factors)
    {
        upolynomial_ZZ cyc24({{24, ZZ(1)}, {0, ZZ(-1)}});
        BENCH("factorize  x^24-1 (cyclotomic)", 3, {
            volatile auto r = factorize(cyc24); (void)r;
        });
    }

    // Swinnerton-Dyer S3 (deg 8, minpoly of sqrt2+sqrt3+sqrt5)
    {
        // S3 = x^8 - 40x^6 + 352x^4 - 960x^2 + 576
        upolynomial_ZZ sd3({{8,ZZ(1)},{6,ZZ(-40)},{4,ZZ(352)},{2,ZZ(-960)},{0,ZZ(576)}});
        BENCH("factorize  Swinnerton-Dyer S3", 5, {
            volatile auto r = factorize(sd3); (void)r;
        });
    }

    // -- squarefree factorization --
    // Build polynomial with repeated factors: f1^1 * f2^2 * f3^3
    auto sf1 = random_polynomial<ZZ>({x}, 7, 6, {-10, 10});
    auto sf2 = random_polynomial<ZZ>({x}, 7, 6, {-10, 10});
    auto sf3 = random_polynomial<ZZ>({x}, 7, 6, {-10, 10});
    auto sqf_s = sf1 * pow(sf2, 2);                      // ~deg21
    auto sqf_m = sf1 * pow(sf2, 2) * pow(sf3, 3);        // ~deg42

    BENCH("squarefree  ~deg21", 10, {
        volatile auto r = squarefreefactorize(sqf_s); (void)r;
    });
    BENCH("squarefree  ~deg42", 5, {
        volatile auto r = squarefreefactorize(sqf_m); (void)r;
    });

    // ================================================================
    // Multivariate polynomial_QQ (2 variables)
    // ================================================================
    bench_header("Multivariate polynomial_QQ (2 variables)");

    auto qq_a_s = random_polynomial_QQ({x, y}, 15, 12, {-30, 30}, 10);
    auto qq_b_s = random_polynomial_QQ({x, y}, 15, 12, {-30, 30}, 10);
    auto qq_a_m = random_polynomial_QQ({x, y}, 30, 25, {-30, 30}, 10);
    auto qq_b_m = random_polynomial_QQ({x, y}, 30, 25, {-30, 30}, 10);

    BENCH("add  deg15", 100, {
        volatile auto r = qq_a_s + qq_b_s; (void)r;
    });
    BENCH("add  deg30", 100, {
        volatile auto r = qq_a_m + qq_b_m; (void)r;
    });

    auto qq_mul_a_s = random_polynomial_QQ({x, y}, 10, 8, {-20, 20}, 10);
    auto qq_mul_b_s = random_polynomial_QQ({x, y}, 10, 8, {-20, 20}, 10);
    auto qq_mul_a_m = random_polynomial_QQ({x, y}, 20, 15, {-20, 20}, 10);
    auto qq_mul_b_m = random_polynomial_QQ({x, y}, 20, 15, {-20, 20}, 10);

    BENCH("mul  deg10", 50, {
        volatile auto r = qq_mul_a_s * qq_mul_b_s; (void)r;
    });
    BENCH("mul  deg20", 10, {
        volatile auto r = qq_mul_a_m * qq_mul_b_m; (void)r;
    });

    auto qq_pow_base = random_polynomial_QQ({x, y}, 3, 4, {-5, 5}, 5);
    BENCH("pow  deg3^5", 20, {
        volatile auto r = pow(qq_pow_base, 5); (void)r;
    });
    BENCH("pow  deg3^8", 5, {
        volatile auto r = pow(qq_pow_base, 8); (void)r;
    });

    // ================================================================
    // Factorize stress tests (many factors)
    // ================================================================
    bench_header("Factorize stress tests");

    // Stress 1: univariate 70 linear factors — (x-1)(x-2)...(x-70)
    {
        upolynomial_ZZ stress_uni({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 70; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            stress_uni = stress_uni * lin;
        }
        BENCH("factorize  uni 70 factors (deg70)", 1, {
            volatile auto r = factorize(stress_uni); (void)r;
        });
    }

    // Stress 2: bivariate 70 linear factors — Π_{i=1}^{35} (x+iy)(x-iy)
    {
        using PolyLex = polynomial_<ZZ, lex>;
        auto to_lex = [](const polynomial_ZZ& p) {
            PolyLex r; poly_convert(p, r); return r;
        };
        auto stress_bi = to_lex(
            (pow(x,1) + pow(y,1)) * (pow(x,1) - pow(y,1)));
        stress_bi.normalization();
        for (int i = 2; i <= 35; ++i) {
            auto f1 = to_lex(pow(x,1) + ZZ(i) * pow(y,1));
            auto f2 = to_lex(pow(x,1) - ZZ(i) * pow(y,1));
            stress_bi = stress_bi * f1;
            stress_bi.normalization();
            stress_bi = stress_bi * f2;
            stress_bi.normalization();
        }
        BENCH("factorize  bivar 70 factors (deg70)", 1, {
            volatile auto r = factorize(stress_bi); (void)r;
        });
    }

    // Stress 3: trivariate 60 linear factors — Π(x±iy, x±iz, y±iz), i=1..10
    {
        variable z("z");
        using PolyLex = polynomial_<ZZ, lex>;
        auto to_lex = [](const polynomial_ZZ& p) {
            PolyLex r; poly_convert(p, r); return r;
        };
        PolyLex stress_tri;
        {
            basic_monomial<lex> m0;
            stress_tri.push_back(std::make_pair(m0, ZZ(1)));
        }
        for (int i = 1; i <= 10; ++i) {
            auto f1 = to_lex(pow(x,1) + ZZ(i) * pow(y,1));
            auto f2 = to_lex(pow(x,1) - ZZ(i) * pow(y,1));
            stress_tri = stress_tri * f1 * f2;
            stress_tri.normalization();
        }
        for (int i = 1; i <= 10; ++i) {
            auto f1 = to_lex(pow(x,1) + ZZ(i) * pow(z,1));
            auto f2 = to_lex(pow(x,1) - ZZ(i) * pow(z,1));
            stress_tri = stress_tri * f1 * f2;
            stress_tri.normalization();
        }
        for (int i = 1; i <= 10; ++i) {
            auto f1 = to_lex(pow(y,1) + ZZ(i) * pow(z,1));
            auto f2 = to_lex(pow(y,1) - ZZ(i) * pow(z,1));
            stress_tri = stress_tri * f1 * f2;
            stress_tri.normalization();
        }
        BENCH("factorize  trivar 60 factors (deg60)", 1, {
            volatile auto r = factorize(stress_tri); (void)r;
        });
    }

    // Stress 4: disjoint pairs — (x1+x2)(x1-x2)*...*(x9+x10)(x9-x10)
    {
        using PolyLex = polynomial_<ZZ, lex>;
        auto to_lex = [](const polynomial_ZZ& p) {
            PolyLex r; poly_convert(p, r); return r;
        };
        std::vector<variable> xvars;
        for (int i = 1; i <= 10; ++i)
            xvars.push_back(variable("x" + std::to_string(i)));

        PolyLex stress_dis;
        {
            basic_monomial<lex> m0;
            stress_dis.push_back(std::make_pair(m0, ZZ(1)));
        }
        for (int i = 0; i < 5; ++i) {
            auto f1 = to_lex(pow(xvars[2*i],1) + pow(xvars[2*i+1],1));
            auto f2 = to_lex(pow(xvars[2*i],1) - pow(xvars[2*i+1],1));
            stress_dis = stress_dis * f1 * f2;
            stress_dis.normalization();
        }
        BENCH("factorize  10var disjoint 10 factors", 1, {
            volatile auto r = factorize(stress_dis); (void)r;
        });
    }

    // ================================================================
    // Real root isolation (uspensky)
    // ================================================================
    bench_header("Real root isolation (uspensky)");

    // Build polynomials with known real roots: (t-1)(t-2)...(t-k)
    {
        // deg 10: (t-1)(t-2)...(t-10)
        upolynomial_ZZ usp_10({{1, ZZ(1)}, {0, ZZ(-1)}});  // t-1
        for (int k = 2; k <= 10; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            usp_10 = usp_10 * lin;
        }
        BENCH("uspensky  deg10 (10 roots)", 20, {
            volatile auto r = uspensky(usp_10); (void)r;
        });

        // deg 20: (t-1)(t-2)...(t-20)
        upolynomial_ZZ usp_20({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 20; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            usp_20 = usp_20 * lin;
        }
        BENCH("uspensky  deg20 (20 roots)", 10, {
            volatile auto r = uspensky(usp_20); (void)r;
        });

        // deg 30: (t-1)(t-2)...(t-30)
        upolynomial_ZZ usp_30({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 30; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            usp_30 = usp_30 * lin;
        }
        BENCH("uspensky  deg30 (30 roots)", 5, {
            volatile auto r = uspensky(usp_30); (void)r;
        });
    }

    // ================================================================
    // Memory Profile: per-polynomial representation size
    // ================================================================
    bench_header("Memory Profile: per-polynomial representation size");
    mem_profile_header();

    const int MEM_N = 500;

    // -- polynomial_ZZ (multivariate) --
    {
        auto p = random_polynomial<ZZ>({x, y}, 20, 15, {-50, 50});
        double kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_ZZ(p); });
        mem_profile_row("polynomial_ZZ deg20 len15", p.size(), kb);
    }
    {
        auto p = random_polynomial<ZZ>({x, y}, 40, 30, {-50, 50});
        double kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_ZZ(p); });
        mem_profile_row("polynomial_ZZ deg40 len30", p.size(), kb);
    }
    {
        auto p = random_polynomial<ZZ>({x, y}, 80, 60, {-50, 50});
        double kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_ZZ(p); });
        mem_profile_row("polynomial_ZZ deg80 len60", p.size(), kb);
    }

    // -- upolynomial_ZZ (univariate) --
    {
        auto p = random_upolynomial<ZZ>(50, 40, {-100, 100});
        double kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
        mem_profile_row("upolynomial_ZZ deg50", p.size(), kb);
    }
    {
        auto p = random_upolynomial<ZZ>(200, 150, {-100, 100});
        double kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
        mem_profile_row("upolynomial_ZZ deg200", p.size(), kb);
    }
    {
        auto p = random_upolynomial<ZZ>(500, 400, {-100, 100});
        double kb = measure_per_object_kb(MEM_N, [&]() { return upolynomial_ZZ(p); });
        mem_profile_row("upolynomial_ZZ deg500", p.size(), kb);
    }

    // -- polynomial_QQ (multivariate) --
    {
        auto p = random_polynomial_QQ({x, y}, 15, 12, {-30, 30}, 10);
        double kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_QQ(p); });
        mem_profile_row("polynomial_QQ deg15", p.size(), kb);
    }
    {
        auto p = random_polynomial_QQ({x, y}, 30, 25, {-30, 30}, 10);
        double kb = measure_per_object_kb(MEM_N, [&]() { return polynomial_QQ(p); });
        mem_profile_row("polynomial_QQ deg30", p.size(), kb);
    }

    print_process_memory();
    std::cout << "\nDone. Total benchmarks: " << _bench_results.size() << std::endl;
    return 0;
}
