/**
 * @file test_open_cad.cc
 * @brief Test open_cad 公开接口
 *
 * 测试策略：只验证 CAD 树的稳定属性（type, parent, sample_point, children），
 * 不验证投影多项式的具体值和 cad_root 中的 poly_idx/root_idx，
 * 因为这些会随分解策略（无平方 vs 因式分解）改变，但不影响 CAD 的正确性。
 * 本测试选取的多项式恰好只有有理数根，realroot 给出精确根，
 * 因此不论分解策略如何，树结构和样本点都是确定的。
 *
 * 全部通过一般序接口 open_cad 调用，同时覆盖了内部 __open_cad 的逻辑。
 */
#include "clpoly_test.hh"
#include <clpoly/cad.hh>
#include <iostream>

using namespace clpoly;

// 验证 Sector 节点的稳定属性（不依赖投影多项式的具体分解方式）
// 只比较 type, parent, sample_point, child_begin, child_count
#define ASSERT_SECTOR(nd, _parent, _sample) \
    CLPOLY_ASSERT_TRUE((nd).is_sector()); \
    CLPOLY_ASSERT_EQ((nd).parent(), size_t(_parent)); \
    CLPOLY_ASSERT_EQ((nd).sample_point(), QQ _sample)

#define ASSERT_SECTOR_CHILDREN(nd, _parent, _sample, _cb, _cc) \
    ASSERT_SECTOR(nd, _parent, _sample); \
    CLPOLY_ASSERT_EQ((nd).child_begin(), size_t(_cb)); \
    CLPOLY_ASSERT_EQ((nd).child_count(), size_t(_cc))

int main() {
    CLPOLY_TEST("open_cad univariate");
    {
        // f = x^2 - 1，实根 x = -1, 1
        // 预期：1 层，3 个 Sector
        //   [0] (-∞, -1)  sample=-2
        //   [1] (-1,  1)  sample=0
        //   [2] ( 1, +∞)  sample=2
        variable x("x");
        using poly_t = polynomial_<ZZ, grlex>;
        poly_t f = poly_t(x) * poly_t(x) - poly_t(1);
        std::vector<poly_t> polys = {f};
        std::vector<variable> vars = {x};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(1));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(3));

        ASSERT_SECTOR(tree.node(0, 0), SIZE_MAX, (-2));
        ASSERT_SECTOR(tree.node(0, 1), SIZE_MAX, (0));
        ASSERT_SECTOR(tree.node(0, 2), SIZE_MAX, (2));
    }

    CLPOLY_TEST("open_cad no real roots");
    {
        // f = x^2 + 1，无实根
        // 预期：1 层，1 个 Sector (-∞, +∞)  sample=0
        variable x("x");
        using poly_t = polynomial_<ZZ, grlex>;
        poly_t f = poly_t(x) * poly_t(x) + poly_t(1);
        std::vector<poly_t> polys = {f};
        std::vector<variable> vars = {x};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(1));
        ASSERT_SECTOR(tree.node(0, 0), SIZE_MAX, (0));
    }

    CLPOLY_TEST("open_cad multiple root");
    {
        // f = (x - 1)^2，重根 x = 1（重数 2）
        // realroot 返回 1 个不同的实根，预期 2 个 Sector
        //   [0] (-∞, 1)  sample=0
        //   [1] (1, +∞)  sample=2
        variable x("x");
        using poly_t = polynomial_<ZZ, grlex>;
        poly_t f = (poly_t(x) - poly_t(1)) * (poly_t(x) - poly_t(1));
        std::vector<poly_t> polys = {f};
        std::vector<variable> vars = {x};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(1));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(2));
        ASSERT_SECTOR(tree.node(0, 0), SIZE_MAX, (0));
        ASSERT_SECTOR(tree.node(0, 1), SIZE_MAX, (2));
    }

    CLPOLY_TEST("open_cad empty polys");
    {
        // 空多项式列表，预期 1 层，1 个 Sector (-∞, +∞)
        variable x("x");
        using poly_t = polynomial_<ZZ, grlex>;
        std::vector<poly_t> polys = {};
        std::vector<variable> vars = {x};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(1));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(1));
        ASSERT_SECTOR(tree.node(0, 0), SIZE_MAX, (0));
    }

    CLPOLY_TEST("open_cad zero polynomial");
    {
        // 零多项式应被 __conts_prims_polys_var 中的 is_number 过滤掉
        // 等价于空多项式列表，预期 1 层，1 个 Sector (-∞, +∞)
        variable x("x");
        using poly_t = polynomial_<ZZ, grlex>;
        poly_t zero;  // 零多项式
        std::vector<poly_t> polys = {zero};
        std::vector<variable> vars = {x};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(1));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(1));
        ASSERT_SECTOR(tree.node(0, 0), SIZE_MAX, (0));
    }

    CLPOLY_TEST("open_cad bivariate");
    {
        // f = x^2 + y^2 - 1（单位圆），vars = {x, y}（lex 序 x < y）
        // 提升序：先 y 后 x
        //
        // Level 0 (y)：投影消去 x 后得到 y^2 - 1 的因子，实根 y = -1, 1
        //   [0] (-∞, -1)  sample=-2  → 圆外，x 无实根 → 1 个子节点
        //   [1] (-1,  1)  sample=0   → 圆内，x^2 = 1-y^2 有 2 根 → 3 个子节点
        //   [2] ( 1, +∞)  sample=2   → 圆外，x 无实根 → 1 个子节点
        //
        // Level 1 (x)：
        //   [0] parent=0  (-∞, +∞)           sample=0   (y=-2, x^2+3=0 无实根)
        //   [1] parent=1  (-∞, -1)           sample=-2  (y=0, x^2-1=0 根为 ±1)
        //   [2] parent=1  (-1,  1)           sample=0
        //   [3] parent=1  ( 1, +∞)           sample=2
        //   [4] parent=2  (-∞, +∞)           sample=0   (y=2, x^2+3=0 无实根)
        variable x("x"), y("y");
        using poly_t = polynomial_<ZZ, grlex>;
        poly_t f = poly_t(x) * poly_t(x) + poly_t(y) * poly_t(y) - poly_t(1);
        std::vector<poly_t> polys = {f};
        std::vector<variable> vars = {x, y};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(2));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(3));
        CLPOLY_ASSERT_EQ(tree.level_size(1), size_t(5));

        // Level 0 (y)
        ASSERT_SECTOR_CHILDREN(tree.node(0, 0), SIZE_MAX, (-2), 0, 1);  // (-∞, -1) 圆外
        ASSERT_SECTOR_CHILDREN(tree.node(0, 1), SIZE_MAX, (0),  1, 3);  // (-1,  1) 圆内
        ASSERT_SECTOR_CHILDREN(tree.node(0, 2), SIZE_MAX, (2),  4, 1);  // ( 1, +∞) 圆外

        // Level 1 (x)
        ASSERT_SECTOR(tree.node(1, 0), 0, (0));    // y=-2 → x^2+3=0 无根，(-∞,+∞)
        ASSERT_SECTOR(tree.node(1, 1), 1, (-2));    // y=0  → x^2-1=0 根±1，(-∞,-1)
        ASSERT_SECTOR(tree.node(1, 2), 1, (0));     // y=0  → x^2-1=0 根±1，(-1, 1)
        ASSERT_SECTOR(tree.node(1, 3), 1, (2));     // y=0  → x^2-1=0 根±1，( 1,+∞)
        ASSERT_SECTOR(tree.node(1, 4), 2, (0));     // y=2  → x^2+3=0 无根，(-∞,+∞)

        for (size_t i = 0; i < tree.level_size(1); ++i)
        {
            auto sp = tree.get_sample_point(1, i);
            CLPOLY_ASSERT_EQ(sp.size(), size_t(2));
        }
    }

    CLPOLY_TEST("open_cad trivariate");
    {
        // f = x^2 + y^2 + z^2 - 1（单位球），vars = {x, y, z}
        // 提升序：先 z 后 y 后 x
        //
        // Level 0 (z)：投影消去 x,y 后得到 z^2 - 1 的因子，实根 z = -1, 1
        //   [0] (-∞, -1) sample=-2  → 球外 → 1 子
        //   [1] (-1,  1) sample=0   → 球内截面是圆 → 3 子
        //   [2] ( 1, +∞) sample=2   → 球外 → 1 子
        //
        // Level 1 (y)：代入 z 的样本点后
        //   [0] parent=0  (-∞,+∞)    sample=0   (z=-2, y^2+3=0 无根 → 1 子)
        //   [1] parent=1  (-∞,-1)    sample=-2  (z=0,  y^2-1=0 根±1, 左区间 → 1 子)
        //   [2] parent=1  (-1, 1)    sample=0   (z=0,  y^2-1=0 根±1, 中区间 → 3 子)
        //   [3] parent=1  ( 1,+∞)    sample=2   (z=0,  y^2-1=0 根±1, 右区间 → 1 子)
        //   [4] parent=2  (-∞,+∞)    sample=0   (z=2,  y^2+3=0 无根 → 1 子)
        //
        // Level 2 (x)：代入 z,y 的样本点后
        //   [0] parent=0  (-∞,+∞)  sample=0   (z=-2,y=0 → x^2+3=0 无根)
        //   [1] parent=1  (-∞,+∞)  sample=0   (z=0,y=-2 → x^2+3=0 无根)
        //   [2] parent=2  (-∞,-1)  sample=-2  (z=0,y=0  → x^2-1=0 根±1)
        //   [3] parent=2  (-1, 1)  sample=0
        //   [4] parent=2  ( 1,+∞)  sample=2
        //   [5] parent=3  (-∞,+∞)  sample=0   (z=0,y=2  → x^2+3=0 无根)
        //   [6] parent=4  (-∞,+∞)  sample=0   (z=2,y=0  → x^2+3=0 无根)
        variable x("x"), y("y"), z("z");
        using poly_t = polynomial_<ZZ, grlex>;
        poly_t f = poly_t(x) * poly_t(x) + poly_t(y) * poly_t(y)
                   + poly_t(z) * poly_t(z) - poly_t(1);
        std::vector<poly_t> polys = {f};
        std::vector<variable> vars = {x, y, z};

        auto tree = open_cad(polys, vars);

        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(3));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(3));
        CLPOLY_ASSERT_EQ(tree.level_size(1), size_t(5));
        CLPOLY_ASSERT_EQ(tree.level_size(2), size_t(7));

        // Level 0 (z)
        ASSERT_SECTOR_CHILDREN(tree.node(0, 0), SIZE_MAX, (-2), 0, 1);  // (-∞,-1) 球外
        ASSERT_SECTOR_CHILDREN(tree.node(0, 1), SIZE_MAX, (0),  1, 3);  // (-1, 1) 球内截面
        ASSERT_SECTOR_CHILDREN(tree.node(0, 2), SIZE_MAX, (2),  4, 1);  // ( 1,+∞) 球外

        // Level 1 (y)
        ASSERT_SECTOR_CHILDREN(tree.node(1, 0), 0, (0),  0, 1);  // z=-2 → y^2+3=0 无根
        ASSERT_SECTOR_CHILDREN(tree.node(1, 1), 1, (-2), 1, 1);  // z=0  → y^2-1=0 根±1，(-∞,-1)
        ASSERT_SECTOR_CHILDREN(tree.node(1, 2), 1, (0),  2, 3);  // z=0  → y^2-1=0 根±1，(-1, 1)
        ASSERT_SECTOR_CHILDREN(tree.node(1, 3), 1, (2),  5, 1);  // z=0  → y^2-1=0 根±1，( 1,+∞)
        ASSERT_SECTOR_CHILDREN(tree.node(1, 4), 2, (0),  6, 1);  // z=2  → y^2+3=0 无根

        // Level 2 (x)
        ASSERT_SECTOR(tree.node(2, 0), 0, (0));    // z=-2,y=0 → x^2+3=0 无根
        ASSERT_SECTOR(tree.node(2, 1), 1, (0));    // z=0,y=-2 → x^2+3=0 无根
        ASSERT_SECTOR(tree.node(2, 2), 2, (-2));    // z=0,y=0  → x^2-1=0 根±1，(-∞,-1)
        ASSERT_SECTOR(tree.node(2, 3), 2, (0));     // z=0,y=0  → x^2-1=0 根±1，(-1, 1)
        ASSERT_SECTOR(tree.node(2, 4), 2, (2));     // z=0,y=0  → x^2-1=0 根±1，( 1,+∞)
        ASSERT_SECTOR(tree.node(2, 5), 3, (0));    // z=0,y=2  → x^2+3=0 无根
        ASSERT_SECTOR(tree.node(2, 6), 4, (0));    // z=2,y=0  → x^2+3=0 无根

        for (size_t i = 0; i < tree.level_size(2); ++i)
        {
            auto sp = tree.get_sample_point(2, i);
            CLPOLY_ASSERT_EQ(sp.size(), size_t(3));
        }
    }

    return clpoly_test::test_summary();
}
