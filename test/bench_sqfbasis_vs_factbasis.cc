/**
 * @file bench_sqfbasis_vs_factbasis.cc
 * @brief 对比 CAD 投影中 squarefreebasis vs factorbasis 的性能
 *
 * 通过 basis_computation_method 开关选择不同的基生成方法
 * 用法: make bench-sqfbasis-vs-factbasis
 */
#include <clpoly/clpoly.hh>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace clpoly;

// 测试用例数据结构
struct TestCase {
    const char* name;
    std::vector<polynomial_ZZ> polys;
    std::vector<variable> vars;
    int proj_repeats;  // project_full 重复次数
    int cad_repeats;   // open_cad 重复次数 (0 = 跳过)
};

struct BenchResult {
    double sqf_proj_ms;
    double fac_proj_ms;
    double sqf_cad_ms;   // -1 = 跳过
    double fac_cad_ms;   // -1 = 跳过
};

// 逐迭代计时，取平均
static BenchResult bench_one(const TestCase& tc)
{
    // lex 转换（md 必须在整个函数存活）
    lex_<custom_var_order> md(tc.vars);
    polynomial_<ZZ, lex_<custom_var_order>> p(&md);
    std::vector<polynomial_<ZZ, lex_<custom_var_order>>> polys_lex;
    polys_lex.reserve(tc.polys.size());
    for (const auto& i : tc.polys) {
        poly_convert(i, p);
        polys_lex.push_back(std::move(p));
    }

    // 预热
    { __project_full(polys_lex, tc.vars, projection_method::LAZARD, basis_computation_method::SQUAREFREE); }
    { __project_full(polys_lex, tc.vars, projection_method::LAZARD, basis_computation_method::FACTOR); }

    // project_full: squarefreebasis
    double sqf_proj = 0;
    for (int r = 0; r < tc.proj_repeats; ++r) {
        auto t0 = std::chrono::high_resolution_clock::now();
        __project_full(polys_lex, tc.vars, projection_method::LAZARD, basis_computation_method::SQUAREFREE);
        auto t1 = std::chrono::high_resolution_clock::now();
        sqf_proj += std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
    sqf_proj /= tc.proj_repeats;

    // project_full: factorbasis
    double fac_proj = 0;
    for (int r = 0; r < tc.proj_repeats; ++r) {
        auto t0 = std::chrono::high_resolution_clock::now();
        __project_full(polys_lex, tc.vars, projection_method::LAZARD, basis_computation_method::FACTOR);
        auto t1 = std::chrono::high_resolution_clock::now();
        fac_proj += std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
    fac_proj /= tc.proj_repeats;

    // open_cad（使用一般序 polys，内部自行转 lex）
    double sqf_cad = -1, fac_cad = -1;
    if (tc.cad_repeats > 0) {
        // 预热
        { open_cad(tc.polys, tc.vars, projection_method::LAZARD, basis_computation_method::SQUAREFREE); }
        { open_cad(tc.polys, tc.vars, projection_method::LAZARD, basis_computation_method::FACTOR); }

        sqf_cad = 0;
        for (int r = 0; r < tc.cad_repeats; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            open_cad(tc.polys, tc.vars, projection_method::LAZARD, basis_computation_method::SQUAREFREE);
            auto t1 = std::chrono::high_resolution_clock::now();
            sqf_cad += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
        sqf_cad /= tc.cad_repeats;

        fac_cad = 0;
        for (int r = 0; r < tc.cad_repeats; ++r) {
            auto t0 = std::chrono::high_resolution_clock::now();
            open_cad(tc.polys, tc.vars, projection_method::LAZARD, basis_computation_method::FACTOR);
            auto t1 = std::chrono::high_resolution_clock::now();
            fac_cad += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }
        fac_cad /= tc.cad_repeats;
    }

    return {sqf_proj, fac_proj, sqf_cad, fac_cad};
}

int main()
{
    variable x("x"), y("y"), z("z");

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "=== CAD: squarefreebasis vs factorbasis ===\n";
#ifdef NDEBUG
    std::cout << "Build: Release\n\n";
#else
    std::cout << "Build: Debug (timings not representative)\n\n";
#endif

    std::cout << std::setw(28) << std::left << "用例"
              << std::setw(22) << std::right << "-- project_full --"
              << "  "
              << std::setw(26) << "--- open_cad ---"
              << "\n";
    std::cout << std::setw(28) << std::left << ""
              << std::setw(9) << std::right << "sqf"
              << std::setw(9) << "fac"
              << std::setw(8) << "ratio"
              << std::setw(9) << "sqf"
              << std::setw(9) << "fac"
              << std::setw(8) << "ratio"
              << "\n";
    std::cout << std::string(86, '-') << "\n";

    auto fmt = [](double ms) -> std::string {
        if (ms < 0) return "---";
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << ms;
        return oss.str();
    };

    auto print_row = [&](const TestCase& tc, const BenchResult& res) {
        auto ratio_str = [](double a, double b) -> std::string {
            if (a < 0 || b < 0) return "---";
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << a / std::max(b, 0.001) << "x";
            return oss.str();
        };
        std::cout << std::setw(28) << std::left << tc.name
                  << std::setw(9) << std::right << fmt(res.sqf_proj_ms)
                  << std::setw(9) << fmt(res.fac_proj_ms)
                  << std::setw(8) << ratio_str(res.sqf_proj_ms, res.fac_proj_ms)
                  << std::setw(9) << fmt(res.sqf_cad_ms)
                  << std::setw(9) << fmt(res.fac_cad_ms)
                  << std::setw(8) << ratio_str(res.sqf_cad_ms, res.fac_cad_ms)
                  << "\n";
    };

    // ---- 经典用例 ----
    std::vector<TestCase> classics = {
        {"circle x^2+y^2-1",                                   
            {pow(x,2)+pow(y,2)-1}, {x,y}, 500, 200},
        {"ellipse x^2+2y^2-1",                                 
            {pow(x,2)+2*pow(y,2)-1}, {x,y}, 500, 200},
        {"{x^2+y^2-1, x-y}",                                   
            {pow(x,2)+pow(y,2)-1, x-y}, {x,y}, 200, 100},
        {"shared factors 2var",                                 
            {(x-1)*(x+1)*(y-1), (x-1)*(y+1)}, {x,y}, 500, 200},
        {"x^4-y^4",                                             
            {pow(x,4)-pow(y,4)}, {x,y}, 500, 200},
        {"sphere x^2+y^2+z^2-1",                               
            {pow(x,2)+pow(y,2)+pow(z,2)-1}, {x,y,z}, 200, 100},
        {"{sphere, plane}",                                     
            {pow(x,2)+pow(y,2)+pow(z,2)-1, x+y+z}, {x,y,z}, 50, 20},
        {"shared factors 3var",                                 
            {(x-1)*(y-z), (x-1)*(y+z), (y-z)*(x+1)}, {x,y,z}, 200, 100},
        {"{x^3+y^3+z^3-1, xyz-1}",                             
            {pow(x,3)+pow(y,3)+pow(z,3)-1, x*y*z-1}, {x,y,z}, 50, 10},
    };

    for (auto& tc : classics) {
        auto res = bench_one(tc);
        print_row(tc, res);
    }

    std::cout << std::string(86, '-') << "\n";

    // ---- 随机二变量 ----
    std::vector<TestCase> randoms = {
        {"rand 2var 1p deg3",                                   
            {random_polynomial<ZZ>({x,y}, 3, 4, {-5,5})},
            {x,y}, 500, 200},
        {"rand 2var 2p deg4",                                   
            {random_polynomial<ZZ>({x,y}, 4, 5, {-5,5}),
             random_polynomial<ZZ>({x,y}, 4, 5, {-5,5})},
            {x,y}, 100, 20},
        {"rand 2var 3p deg5",                                   
            {random_polynomial<ZZ>({x,y}, 5, 6, {-5,5}),
             random_polynomial<ZZ>({x,y}, 5, 6, {-5,5}),
             random_polynomial<ZZ>({x,y}, 5, 6, {-5,5})},
            {x,y}, 20, 5},
        {"rand 2var 3p deg8",                                   
            {random_polynomial<ZZ>({x,y}, 8, 10, {-5,5}),
             random_polynomial<ZZ>({x,y}, 8, 10, {-5,5}),
             random_polynomial<ZZ>({x,y}, 8, 10, {-5,5})},
            {x,y}, 5, 3},
        {"rand 2var 5p deg6",                                   
            {random_polynomial<ZZ>({x,y}, 6, 8, {-5,5}),
             random_polynomial<ZZ>({x,y}, 6, 8, {-5,5}),
             random_polynomial<ZZ>({x,y}, 6, 8, {-5,5}),
             random_polynomial<ZZ>({x,y}, 6, 8, {-5,5}),
             random_polynomial<ZZ>({x,y}, 6, 8, {-5,5})},
            {x,y}, 5, 3},
    };

    for (auto& tc : randoms) {
        auto res = bench_one(tc);
        print_row(tc, res);
    }

    std::cout << std::string(86, '-') << "\n";

    // ---- 随机三变量 ----
    std::vector<TestCase> randoms_3var = {
        {"rand 3var 1p deg3",                                   
            {random_polynomial<ZZ>({x,y,z}, 3, 5, {-5,5})},
            {x,y,z}, 200, 50},
        {"rand 3var 2p deg3",                                   
            {random_polynomial<ZZ>({x,y,z}, 3, 5, {-5,5}),
             random_polynomial<ZZ>({x,y,z}, 3, 5, {-5,5})},
            {x,y,z}, 10, 3},
        {"rand 3var 2p deg4",                                   
            {random_polynomial<ZZ>({x,y,z}, 4, 6, {-5,5}),
             random_polynomial<ZZ>({x,y,z}, 4, 6, {-5,5})},
            {x,y,z}, 5, 3},
        {"rand 3var 3p deg3",                                   
            {random_polynomial<ZZ>({x,y,z}, 3, 5, {-5,5}),
             random_polynomial<ZZ>({x,y,z}, 3, 5, {-5,5}),
             random_polynomial<ZZ>({x,y,z}, 3, 5, {-5,5})},
            {x,y,z}, 5, 3},
    };

    for (auto& tc : randoms_3var) {
        auto res = bench_one(tc);
        print_row(tc, res);
    }

    std::cout << "\n";
    return 0;
}
