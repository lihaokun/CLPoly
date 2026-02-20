#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>
#include <chrono>
#include <vector>

using namespace clpoly;
using namespace std;

struct Timer {
    chrono::high_resolution_clock::time_point start;
    Timer() : start(chrono::high_resolution_clock::now()) {}
    double elapsed_ms() {
        auto now = chrono::high_resolution_clock::now();
        return chrono::duration<double, milli>(now - start).count();
    }
};

int main()
{
    // ======================================================
    // Test 1: 单变量 70 因子 — 测试 r>64 重组
    // f = (x-1)(x-2)...(x-70)
    // ======================================================
    {
        cout << "=== Test 1: univariate 70 linear factors ===" << endl;

        upolynomial_<ZZ> f({{umonomial(0), ZZ(1)}});
        for (int i = 1; i <= 70; ++i)
        {
            upolynomial_<ZZ> factor({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-i)}});
            f = f * factor;
            f.normalization();
        }
        cout << "  deg(f) = " << get_deg(f) << ", terms = " << f.size() << endl;

        Timer t;
        auto fac = factorize(f);
        double ms = t.elapsed_ms();

        cout << "  factors found: " << fac.factors.size()
             << ", content = " << fac.content
             << ", time = " << ms << " ms" << endl;

        upolynomial_<ZZ> check({{umonomial(0), fac.content}});
        for (auto& [fi, ei] : fac.factors)
            for (uint64_t e = 0; e < ei; ++e)
            {
                check = check * fi;
                check.normalization();
            }
        cout << "  verification: " << (check == f ? "PASS" : "FAIL") << endl;
        cout << endl;
    }

    // ======================================================
    // Test 2: 二变量 70 因子
    // f = Π_{i=1}^{35} (x + i*y)(x - i*y) = Π(x^2 - i^2*y^2)
    // 70 个线性因子, 2 变量, ~36 项
    // ======================================================
    {
        cout << "=== Test 2: bivariate 70 linear factors ===" << endl;

        variable x("x"), y("y");
        polynomial_ZZ product = pow(x, 1) + pow(y, 1);  // (x + y)
        product = product * (pow(x, 1) - pow(y, 1));     // * (x - y)
        product.normalization();

        for (int i = 2; i <= 35; ++i)
        {
            polynomial_ZZ f1 = pow(x, 1) + ZZ(i) * pow(y, 1);
            polynomial_ZZ f2 = pow(x, 1) - ZZ(i) * pow(y, 1);
            product = product * f1;
            product.normalization();
            product = product * f2;
            product.normalization();
        }
        cout << "  terms = " << product.size()
             << ", deg = " << product.front().first.deg() << endl;

        Timer t;
        auto fac = factorize(product);
        double ms = t.elapsed_ms();

        cout << "  factors found: " << fac.factors.size()
             << ", time = " << ms << " ms" << endl;

        // 验证
        polynomial_ZZ check;
        {
            basic_monomial<grlex> m0;
            check.push_back(std::make_pair(m0, fac.content));
        }
        for (auto& [fi, ei] : fac.factors)
            for (uint64_t e = 0; e < ei; ++e)
            {
                check = check * fi;
                check.normalization();
            }
        cout << "  verification: " << (check == product ? "PASS" : "FAIL") << endl;
        cout << endl;
    }

    // ======================================================
    // Test 3: 三变量 多因子
    // f = Π_{i=1}^{10} (x+iy)(x-iy) * Π_{i=1}^{10} (x+iz)(x-iz) * Π_{i=1}^{10} (y+iz)(y-iz)
    // 60 个线性因子, 3 变量
    // ======================================================
    {
        cout << "=== Test 3: trivariate 60 linear factors ===" << endl;

        variable x("x"), y("y"), z("z");

        polynomial_ZZ product;
        {
            basic_monomial<grlex> m0;
            product.push_back(std::make_pair(m0, ZZ(1)));
        }

        int nfactors = 0;
        // (x ± i*y) for i=1..10
        for (int i = 1; i <= 10; ++i)
        {
            polynomial_ZZ f1 = pow(x, 1) + ZZ(i) * pow(y, 1);
            polynomial_ZZ f2 = pow(x, 1) - ZZ(i) * pow(y, 1);
            product = product * f1 * f2;
            product.normalization();
            nfactors += 2;
        }
        cout << "  after x±iy: " << product.size() << " terms, "
             << nfactors << " factors" << endl;

        // (x ± i*z) for i=1..10
        for (int i = 1; i <= 10; ++i)
        {
            polynomial_ZZ f1 = pow(x, 1) + ZZ(i) * pow(z, 1);
            polynomial_ZZ f2 = pow(x, 1) - ZZ(i) * pow(z, 1);
            product = product * f1 * f2;
            product.normalization();
            nfactors += 2;
        }
        cout << "  after x±iz: " << product.size() << " terms, "
             << nfactors << " factors" << endl;

        // (y ± i*z) for i=1..10
        for (int i = 1; i <= 10; ++i)
        {
            polynomial_ZZ f1 = pow(y, 1) + ZZ(i) * pow(z, 1);
            polynomial_ZZ f2 = pow(y, 1) - ZZ(i) * pow(z, 1);
            product = product * f1 * f2;
            product.normalization();
            nfactors += 2;
        }
        cout << "  final: " << product.size() << " terms, "
             << nfactors << " factors, deg = "
             << product.front().first.deg() << endl;

        if (product.size() <= 500000)
        {
            Timer t;
            auto fac = factorize(product);
            double ms = t.elapsed_ms();
            cout << "  factors found: " << fac.factors.size()
                 << ", time = " << ms << " ms" << endl;
        }
        else
        {
            cout << "  [SKIP] too many terms" << endl;
        }
        cout << endl;
    }

    // ======================================================
    // Test 4: 不相交变量对, 控制规模
    // (x1+x2)(x1-x2) * (x3+x4)(x3-x4) * ...
    // 每对贡献 2 因子 + 2 项, 总项 = 2^(对数)
    // 5 对 = 10 因子, 10 变量, 32 项
    // ======================================================
    {
        cout << "=== Test 4: disjoint pairs, 10 vars, 10 factors ===" << endl;

        vector<variable> vars;
        for (int i = 1; i <= 10; ++i)
            vars.push_back(variable("x" + to_string(i)));

        polynomial_ZZ product;
        {
            basic_monomial<grlex> m0;
            product.push_back(std::make_pair(m0, ZZ(1)));
        }

        for (int i = 0; i < 5; ++i)
        {
            polynomial_ZZ f1 = pow(vars[2*i], 1) + pow(vars[2*i+1], 1);
            polynomial_ZZ f2 = pow(vars[2*i], 1) - pow(vars[2*i+1], 1);
            product = product * f1 * f2;
            product.normalization();
        }
        cout << "  terms = " << product.size()
             << ", deg = " << product.front().first.deg() << endl;

        Timer t;
        auto fac = factorize(product);
        double ms = t.elapsed_ms();

        cout << "  factors found: " << fac.factors.size()
             << ", time = " << ms << " ms" << endl;

        // 验证
        polynomial_ZZ check;
        {
            basic_monomial<grlex> m0;
            check.push_back(std::make_pair(m0, fac.content));
        }
        for (auto& [fi, ei] : fac.factors)
            for (uint64_t e = 0; e < ei; ++e)
            {
                check = check * fi;
                check.normalization();
            }
        cout << "  verification: " << (check == product ? "PASS" : "FAIL") << endl;
        cout << endl;
    }

    return 0;
}
