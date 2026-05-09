/**
 * @file test_squarefreebasis_bugfix.cc
 * @brief Issue #14 回归测试 + sqf canonical 形式断言 + sanity 用例
 *
 * 修正方案：docs/fixes/2026-05-09-squarefreebasis-unit-leak.md
 *
 * Bug：clpoly/polynomial_gcd.hh:109-111 squarefreefactorize 未做首项 sign 归一化，
 * 导致 factor 与 polynomial_GCD 输出在 associates 上字面不等，使 squarefreebasis
 * 在拆分时把 unit ±1 残留进 basis。
 *
 * 修复：sqf 的 F_ 在 cont 提取后再做一次"首项 sign 提到 f_cont"，使 factor
 * 永远首项正。对齐 SymPy/Maple/FLINT 约定。
 */
#include <clpoly/clpoly.hh>
#include <iostream>
#include <vector>

using namespace clpoly;

static int passed = 0, total = 0;

#define CHECK(cond, label) do {                                     \
    ++total;                                                        \
    if (cond) { ++passed; std::cout << "  [PASS] " << label << std::endl; } \
    else { std::cout << "  [FAIL] " << label << std::endl; }        \
} while(0)

template <class P>
static bool basis_no_unit(const std::vector<P>& basis)
{
    for (auto& p : basis) if (is_number(p)) return false;
    return true;
}

template <class P>
static bool basis_all_squarefree(const std::vector<P>& basis)
{
    for (auto& p : basis) if (!is_number(p) && !is_squarefree(p)) return false;
    return true;
}

template <class P>
static void show(const std::string& label, const std::vector<P>& basis)
{
    std::cout << "    " << label << " basis: [";
    for (size_t i = 0; i < basis.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << basis[i];
    }
    std::cout << "]" << std::endl;
}

int main()
{
    variable x("x"), y("y");

    // ============================================================
    // 主回归：issue #14 复现用例（关键）
    // ============================================================
    std::cout << "=== Issue #14 回归 ===" << std::endl;
    {
        std::vector<polynomial_<ZZ>> polys = {-x+1, (x-1)*x};
        auto [basis, _] = squarefreebasis(polys);
        show("issue14_ex1", basis);
        CHECK(basis_no_unit(basis), "issue14_ex1: basis 无 unit 元素");
    }
    {
        std::vector<polynomial_<ZZ>> polys = {-x+1, (x-1)*x, -x+2, (x-2)*x};
        auto [basis, _] = squarefreebasis(polys);
        show("issue14_ex2", basis);
        CHECK(basis_no_unit(basis), "issue14_ex2: basis 无 unit 元素");
    }
    {
        std::vector<polynomial_<ZZ>> polys = {-x+1, (x-1)*x, -x+2, (x-2)*x,
                                                -x+3, (x-3)*x};
        auto [basis, _] = squarefreebasis(polys);
        show("issue14_ex3", basis);
        CHECK(basis_no_unit(basis), "issue14_ex3: basis 无 unit 元素");
    }

    // ============================================================
    // sqf canonical 形式断言：非常数 factor 首项系数 > 0
    // ============================================================
    std::cout << "\n=== sqf canonical 形式断言 ===" << std::endl;
    {
        polynomial_<ZZ> f = -x + 1;
        auto sqf = squarefreefactorize(f);
        bool ok = true;
        // 跳过第一项（cont），后续 factor 必须首项正
        for (auto it = sqf.begin() + 1; it != sqf.end(); ++it) {
            if (!is_number(it->first) && it->first.front().second < ZZ(0)) {
                ok = false; break;
            }
        }
        std::cout << "    sqf(-x+1) = ";
        for (auto& [p, e] : sqf) std::cout << "(" << p << "," << e << ") ";
        std::cout << std::endl;
        CHECK(ok, "sqf_neg_lead_canonical: -x+1 的非常数 factor 首项 > 0");
    }
    {
        polynomial_<ZZ> f = -6 * (x-1) * (x-2);
        auto sqf = squarefreefactorize(f);
        bool ok = true;
        for (auto it = sqf.begin() + 1; it != sqf.end(); ++it) {
            if (!is_number(it->first) && it->first.front().second < ZZ(0)) {
                ok = false; break;
            }
        }
        std::cout << "    sqf(-6(x-1)(x-2)) = ";
        for (auto& [p, e] : sqf) std::cout << "(" << p << "," << e << ") ";
        std::cout << std::endl;
        CHECK(ok, "sqf_neg_multivar: -6(x-1)(x-2) 非常数 factor 首项 > 0");
    }
    {
        // 多变量 + 首项负
        polynomial_<ZZ> f = -x*y + 1;
        auto sqf = squarefreefactorize(f);
        bool ok = true;
        for (auto it = sqf.begin() + 1; it != sqf.end(); ++it) {
            if (!is_number(it->first) && it->first.front().second < ZZ(0)) {
                ok = false; break;
            }
        }
        std::cout << "    sqf(-x*y+1) = ";
        for (auto& [p, e] : sqf) std::cout << "(" << p << "," << e << ") ";
        std::cout << std::endl;
        CHECK(ok, "sqf_neg_multivar_xy: -xy+1 非常数 factor 首项 > 0");
    }
    {
        // 不变量：sqf 乘积重建 = 原 F（含 cont * 各 factor^mult）
        polynomial_<ZZ> f = -6 * (x-1) * (x-2);
        auto sqf = squarefreefactorize(f);
        polynomial_<ZZ> reconstructed({{monomial(), ZZ(1)}});
        for (auto& [p, e] : sqf) reconstructed = reconstructed * pow(p, e);
        std::cout << "    rebuild f=" << f << " from sqf, got=" << reconstructed
                  << std::endl;
        CHECK(reconstructed == f, "sqf_product_invariant: ∏ pᵢ^eᵢ = F");
    }

    // ============================================================
    // Sanity（不破坏正向用例）
    // ============================================================
    std::cout << "\n=== Sanity（修复不破坏其他正向用例） ===" << std::endl;
    {
        std::vector<polynomial_<ZZ>> polys = {x+1, (x+1)*x};
        auto [basis, _] = squarefreebasis(polys);
        show("sanity_pos_lead", basis);
        CHECK(basis.size() == 2 && basis_no_unit(basis),
              "sanity_pos_lead: {x+1, (x+1)*x} 无符号歧义场景仍正确");
    }
    {
        std::vector<polynomial_<ZZ>> polys = {(x-1)*(x-1), x*(x-1)};
        auto [basis, idx] = squarefreebasis(polys);
        show("sanity_repeated_factor", basis);
        // 期望：basis 包含 (x-1) 和 x 的 squarefree 表示，无 unit
        // multiplicity 在 idx 中正确（x-1 在 F[0] 出现 2 次）
        bool mult_ok = false;
        for (size_t k = 0; k < basis.size(); ++k) {
            // 找 x-1 对应的 basis 元素，看 idx[k] 是否记录 (0, 2)
            for (auto& [i, e] : idx[k]) {
                if (i == 0 && e == 2) mult_ok = true;
            }
        }
        CHECK(basis_no_unit(basis) && basis_all_squarefree(basis) && mult_ok,
              "sanity_repeated_factor: 重数 2 正确捕获，无 unit");
    }
    {
        std::vector<polynomial_<ZZ>> polys = {polynomial_<ZZ>(x), x+1};
        auto [basis, _] = squarefreebasis(polys);
        show("sanity_coprime", basis);
        CHECK(basis.size() == 2 && basis_no_unit(basis),
              "sanity_coprime: {x, x+1} 互素输入返回原集合");
    }

    // ============================================================
    // Summary
    // ============================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary: " << passed << "/" << total << " passed" << std::endl;
    return passed == total ? 0 : 1;
}
