/**
 * @file test_cad_projector.cc
 * @brief Test projection operators (MCCALLUM and LAZARD)
 *
 * This test verifies that projection computations run without errors.
 */
#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <iostream>

using namespace clpoly;

// 计算squarefreebasis的prod
polynomial_ZZ squarefreebasis_prod(const std::vector<polynomial_ZZ>& F)
{
    auto basis_info = squarefreebasis(F);
    const auto& basis = basis_info.first;  // list of squarefree factors
    polynomial_ZZ prod = 1;
    for (const auto& b : basis) {
        prod = prod * b;
    }
    return prod;
}

// 检查squarefreebasis的prod是否相等
void verify_proj_eq(const std::vector<polynomial_ZZ>& F,const std::vector<polynomial_ZZ>& G)
{
    polynomial_ZZ prodF = pow(squarefreebasis_prod(F),2);
    polynomial_ZZ prodG = pow(squarefreebasis_prod(G),2);
    CLPOLY_ASSERT_EQ(prodF, prodG);
}

int main() {
    // Define variables and lexicographic order (x1 < x2 < x3)
    variable x1("x1"), x2("x2"), x3("x3");
    // Construct polynomials
    polynomial_ZZ p1, p2, p3;
    p1 = x3*x3 + x2*x3 + x1 - 1;
    p2 = x2*x2 - x1*x3 + 2;
    p3 = x3 - x1*x1 - x2;

    std::vector<polynomial_ZZ> polys = {p1, p2, p3};

    // ----- MCCALLUM projection -----
    CLPOLY_TEST("MCCALLUM Projection operators");
    {
        auto mccallum_res = polys;
        std::vector<polynomial_ZZ> expected_res;
        // Project x3
        mccallum_res = project(mccallum_res, x3, projection_method::MCCALLUM);
        expected_res = {x1-1, x2, pow(x2,2)+polynomial_ZZ({{{},2}}), x1, \
pow(x1,2)+x2, -pow(x2,2)+4*x1+polynomial_ZZ({{{},-4}}), \
pow(x2,4)+pow(x2,3)*x1+pow(x1,3)+4*pow(x2,2)+2*x2*x1-pow(x1,2)+polynom\
ial_ZZ({{{},4}}), pow(x1,4)+3*x2*pow(x1,2)+2*pow(x2,2)+x1-1, \
-pow(x1,3)+pow(x2,2)-x2*x1+polynomial_ZZ({{{},2}})};
        verify_proj_eq(mccallum_res,expected_res);

        // Project x2
        mccallum_res = project(mccallum_res, x2, projection_method::MCCALLUM);
        expected_res = {x1-1, x1, polynomial_ZZ(), polynomial_ZZ({{{},2}}), \
polynomial_ZZ({{{},4}}), pow(x1,3)-pow(x1,2)+polynomial_ZZ({{{},4}}), \
pow(x1,4)+x1-1, polynomial_ZZ({{{},3}}), \
pow(x1,3)+polynomial_ZZ({{{},-2}}), polynomial_ZZ({{{},-8}}), \
polynomial_ZZ({{{},16}}), \
27*pow(x1,6)-310*pow(x1,5)+603*pow(x1,4)-288*pow(x1,3)-1024*pow(x1,2)+1\
024*x1+polynomial_ZZ({{{},-256}}), \
pow(x1,4)-8*x1+polynomial_ZZ({{{},8}}), \
4*pow(x1,3)+pow(x1,2)+polynomial_ZZ({{{},-8}}), \
pow(x1,4)+polynomial_ZZ({{{},2}}), 2*x1-1, \
pow(x1,8)+2*pow(x1,5)+8*pow(x1,4)+pow(x1,2)-10*x1+polynomial_\
ZZ({{{},25}}), pow(x1,4)-4*x1+polynomial_ZZ({{{},4}}), \
pow(x1,8)-pow(x1,7)+4*pow(x1,4)-pow(x1,3)-pow(x1,2)+polynomial_\
ZZ({{{},4}}), pow(x1,3)-17*pow(x1,2)+16*x1+polynomial_ZZ({{{},-4}}), \
pow(x1,4)-9*x1+polynomial_ZZ({{{},9}}), \
pow(x1,6)-8*pow(x1,4)+20*pow(x1,2)-16*x1+polynomial_ZZ({{{},4}}), \
pow(x1,8)-2*pow(x1,7)+2*pow(x1,5)+15*pow(x1,4)-9*pow(x1,3)-pow(x1,2)-1\
0*x1+polynomial_ZZ({{{},25}}), \
pow(x1,8)-pow(x1,7)+2*pow(x1,5)-2*pow(x1,4)+10*pow(x1,3)-7*pow(x1,2)-1\
0*x1+polynomial_ZZ({{{},25}})};
        verify_proj_eq(mccallum_res,expected_res);
    }

    // ----- LAZARD projection -----
    CLPOLY_TEST("LAZARD Projection operators");
    {
        auto lazard_res = polys;   // intentionally keeping the typo
        std::vector<polynomial_ZZ> expected_res;
        // Project x3
        lazard_res = project(lazard_res, x3, projection_method::LAZARD);
        expected_res = {x1-1, pow(x2,2)+polynomial_ZZ({{{},2}}), x1, pow(x1,2)+x2, \
-pow(x2,2)+4*x1+polynomial_ZZ({{{},-4}}), \
pow(x2,4)+pow(x2,3)*x1+pow(x1,3)+4*pow(x2,2)+2*x2*x1-pow(x1,2)+polynom\
ial_ZZ({{{},4}}), pow(x1,4)+3*x2*pow(x1,2)+2*pow(x2,2)+x1-1, \
-pow(x1,3)+pow(x2,2)-x2*x1+polynomial_ZZ({{{},2}})};
        verify_proj_eq(lazard_res,expected_res);

        // Project x2
        lazard_res = project(lazard_res, x2, projection_method::LAZARD);
        expected_res = {x1-1, x1, polynomial_ZZ({{{},2}}), polynomial_ZZ({{{},4}}), \
pow(x1,3)-pow(x1,2)+polynomial_ZZ({{{},4}}), pow(x1,4)+x1-1, \
pow(x1,3)+polynomial_ZZ({{{},-2}}), polynomial_ZZ({{{},-8}}), \
polynomial_ZZ({{{},16}}), \
27*pow(x1,6)-310*pow(x1,5)+603*pow(x1,4)-288*pow(x1,3)-1024*pow(x1,2)+1\
024*x1+polynomial_ZZ({{{},-256}}), \
pow(x1,4)-8*x1+polynomial_ZZ({{{},8}}), \
4*pow(x1,3)+pow(x1,2)+polynomial_ZZ({{{},-8}}), \
pow(x1,4)+polynomial_ZZ({{{},2}}), 2*x1-1, \
pow(x1,8)+2*pow(x1,5)+8*pow(x1,4)+pow(x1,2)-10*x1+polynomial_\
ZZ({{{},25}}), pow(x1,4)-4*x1+polynomial_ZZ({{{},4}}), \
pow(x1,8)-pow(x1,7)+4*pow(x1,4)-pow(x1,3)-pow(x1,2)+polynomial_\
ZZ({{{},4}}), pow(x1,3)-17*pow(x1,2)+16*x1+polynomial_ZZ({{{},-4}}), \
pow(x1,4)-9*x1+polynomial_ZZ({{{},9}}), \
pow(x1,6)-8*pow(x1,4)+20*pow(x1,2)-16*x1+polynomial_ZZ({{{},4}}), \
pow(x1,8)-2*pow(x1,7)+2*pow(x1,5)+15*pow(x1,4)-9*pow(x1,3)-pow(x1,2)-1\
0*x1+polynomial_ZZ({{{},25}}), \
pow(x1,8)-pow(x1,7)+2*pow(x1,5)-2*pow(x1,4)+10*pow(x1,3)-7*pow(x1,2)-1\
0*x1+polynomial_ZZ({{{},25}})};
        verify_proj_eq(lazard_res,expected_res);
    }

    return clpoly_test::test_summary();
}
