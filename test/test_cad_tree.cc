/**
 * @file test_cad_tree.cc
 * @brief Test cad_root, cad_node, cad_tree 数据结构
 */
#include "clpoly_test.hh"
#include <clpoly/cad_tree.hh>
#include <iostream>

using namespace clpoly;

int main() {
    // ===== cad_root =====
    CLPOLY_TEST("cad_root creation");
    {
        auto r_inf = cad_root::inf();
        CLPOLY_ASSERT_TRUE(r_inf.isinf());
        CLPOLY_ASSERT_FALSE(r_inf.isneginf());
        CLPOLY_ASSERT_FALSE(r_inf.isfinite());

        auto r_neginf = cad_root::neginf();
        CLPOLY_ASSERT_TRUE(r_neginf.isneginf());
        CLPOLY_ASSERT_FALSE(r_neginf.isinf());
        CLPOLY_ASSERT_FALSE(r_neginf.isfinite());

        auto r = cad_root(2, 3);
        CLPOLY_ASSERT_TRUE(r.isfinite());
        CLPOLY_ASSERT_EQ(r.poly_idx(), size_t(2));
        CLPOLY_ASSERT_EQ(r.root_idx(), size_t(3));

        // std::cout << "  output: " << r_neginf << " " << r << " " << r_inf << std::endl;
    }

    // ===== cad_node =====
    CLPOLY_TEST("cad_node SECTOR");
    {
        auto left = cad_root::neginf();
        auto right = cad_root(0, 0);
        cad_node sector(SIZE_MAX, left, right, QQ(1, 2));

        CLPOLY_ASSERT_TRUE(sector.is_sector());
        CLPOLY_ASSERT_FALSE(sector.is_section());
        CLPOLY_ASSERT_EQ(sector.parent(), SIZE_MAX);
        CLPOLY_ASSERT_EQ(sector.roots().size(), size_t(2));
        CLPOLY_ASSERT_TRUE(sector.roots()[0].isneginf());
        CLPOLY_ASSERT_TRUE(sector.roots()[1].isfinite());
        CLPOLY_ASSERT_EQ(sector.sample_point(), QQ(1, 2));
        CLPOLY_ASSERT_EQ(sector.child_count(), size_t(0));

        // std::cout << "  output: " << sector << std::endl;
    }

    CLPOLY_TEST("cad_node SECTION");
    {
        auto root = cad_root(1, 0);
        cad_node section(5, root);

        CLPOLY_ASSERT_TRUE(section.is_section());
        CLPOLY_ASSERT_FALSE(section.is_sector());
        CLPOLY_ASSERT_EQ(section.parent(), size_t(5));
        CLPOLY_ASSERT_EQ(section.roots().size(), size_t(1));
        CLPOLY_ASSERT_EQ(section.roots()[0].poly_idx(), size_t(1));
        CLPOLY_ASSERT_EQ(section.roots()[0].root_idx(), size_t(0));

        // std::cout << "  output: " << section << std::endl;
    }

    // ===== cad_tree =====
    CLPOLY_TEST("cad_tree construction");
    {
        // 空树
        cad_tree<less> empty_tree;
        CLPOLY_ASSERT_EQ(empty_tree.num_levels(), size_t(0));

        // 带变量的树
        variable x("x"), y("y");
        using poly_type = polynomial_<ZZ, lex>;
        std::vector<variable> vars = {x, y};
        std::vector<std::vector<poly_type>> allprojs(2);

        cad_tree<less> tree(vars, std::move(allprojs));
        CLPOLY_ASSERT_EQ(tree.num_levels(), size_t(2));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(0));
        CLPOLY_ASSERT_EQ(tree.level_size(1), size_t(0));
        // 提升序: level 0 = y (vars[1])，level 1 = x (vars[0])
        CLPOLY_ASSERT_EQ(tree.level_var(0), y);
        CLPOLY_ASSERT_EQ(tree.level_var(1), x);
    }

    CLPOLY_TEST("cad_tree add_sector and parent/children");
    {
        variable x("x"), y("y");
        using poly_type = polynomial_<ZZ, lex>;
        std::vector<variable> vars = {x, y};
        std::vector<std::vector<poly_type>> allprojs(2);

        cad_tree<less> tree(vars, std::move(allprojs));

        // Level 0: 3 个 sector（无 parent）
        auto idx0 = tree.add_sector(0, SIZE_MAX, cad_root::neginf(), cad_root(0, 0), QQ(-2));
        auto idx1 = tree.add_sector(0, SIZE_MAX, cad_root(0, 0), cad_root(0, 1), QQ(0));
        auto idx2 = tree.add_sector(0, SIZE_MAX, cad_root(0, 1), cad_root::inf(), QQ(3));

        CLPOLY_ASSERT_EQ(idx0, size_t(0));
        CLPOLY_ASSERT_EQ(idx1, size_t(1));
        CLPOLY_ASSERT_EQ(idx2, size_t(2));
        CLPOLY_ASSERT_EQ(tree.level_size(0), size_t(3));

        // Level 1: 在 Level 0 第 0 个节点下添加子节点
        auto idx10 = tree.add_sector(1, 0, cad_root::neginf(), cad_root::inf(), QQ(1));
        CLPOLY_ASSERT_EQ(idx10, size_t(0));
        CLPOLY_ASSERT_EQ(tree.level_size(1), size_t(1));

        // 检查 parent/children 关系
        CLPOLY_ASSERT_EQ(tree.node(1, 0).parent(), size_t(0));
        CLPOLY_ASSERT_EQ(tree.node(0, 0).child_count(), size_t(1));
        CLPOLY_ASSERT_EQ(tree.node(0, 0).child_begin(), size_t(0));
        // 其他 level 0 节点无 children
        CLPOLY_ASSERT_EQ(tree.node(0, 1).child_count(), size_t(0));
        CLPOLY_ASSERT_EQ(tree.node(0, 2).child_count(), size_t(0));
    }

    CLPOLY_TEST("cad_tree get_sample_point");
    {
        variable x("x"), y("y"), z("z");
        using poly_type = polynomial_<ZZ, lex>;
        std::vector<variable> vars = {x, y, z};
        std::vector<std::vector<poly_type>> allprojs(3);

        cad_tree<less> tree(vars, std::move(allprojs));

        // Level 0: sample = 5
        tree.add_sector(0, SIZE_MAX, cad_root::neginf(), cad_root::inf(), QQ(5));
        // Level 1: sample = -3
        tree.add_sector(1, 0, cad_root::neginf(), cad_root::inf(), QQ(-3));
        // Level 2: sample = 7
        tree.add_sector(2, 0, cad_root::neginf(), cad_root::inf(), QQ(7));

        auto path = tree.get_sample_point(2, 0);
        CLPOLY_ASSERT_EQ(path.size(), size_t(3));
        CLPOLY_ASSERT_EQ(path[0], QQ(5));
        CLPOLY_ASSERT_EQ(path[1], QQ(-3));
        CLPOLY_ASSERT_EQ(path[2], QQ(7));
    }

    // CLPOLY_TEST("cad_tree print");
    // {
    //     variable x("x"), y("y");
    //     using poly_type = polynomial_<ZZ, lex>;
    //     std::vector<variable> vars = {x, y};
    //     std::vector<std::vector<poly_type>> allprojs(2);

    //     cad_tree<less> tree(vars, std::move(allprojs));

    //     tree.add_sector(0, SIZE_MAX, cad_root::neginf(), cad_root(0, 0), QQ(-1));
    //     tree.add_sector(0, SIZE_MAX, cad_root(0, 0), cad_root::inf(), QQ(2));
    //     tree.add_sector(1, 0, cad_root::neginf(), cad_root::inf(), QQ(0));
    //     tree.add_sector(1, 1, cad_root::neginf(), cad_root::inf(), QQ(0));

    //     std::cout << tree;
    // }

    return clpoly_test::test_summary();
}
