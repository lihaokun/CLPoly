#include <clpoly/clpoly.hh>
#include <clpoly/associatedgraph.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    // ======== graph: basic construction ========
    CLPOLY_TEST("graph_basic");
    {
        graph<variable> G;
        variable x("gx"), y("gy"), z("gz");
        G.add_node(x);
        G.add_node(y);
        G.add_node(z);
        CLPOLY_ASSERT_EQ((int)G.size(), 3);

        G.add_edge(x, y);
        CLPOLY_ASSERT_TRUE(G.is_edge(x, y));
        CLPOLY_ASSERT_TRUE(G.is_edge(y, x));  // undirected
        CLPOLY_ASSERT_FALSE(G.is_edge(x, z));
    }

    // ======== associatedgraph ========
    CLPOLY_TEST("associatedgraph_basic");
    {
        variable x1("ag1"), x2("ag2"), x3("ag3"), x4("ag4");
        // f1 = x1*x4 - x2*x3 involves x1,x2,x3,x4
        polynomial_ZZ f1 = x1*x4 - x2*x3;
        std::vector<polynomial_ZZ> F = {f1};
        auto G = associatedgraph(F);

        // All 4 variables should be nodes
        CLPOLY_ASSERT_EQ((int)G.size(), 4);
        // x1-x2, x1-x3, x1-x4, x2-x3, x2-x4, x3-x4 all connected
        CLPOLY_ASSERT_TRUE(G.is_edge(x1, x2));
        CLPOLY_ASSERT_TRUE(G.is_edge(x1, x4));
        CLPOLY_ASSERT_TRUE(G.is_edge(x2, x3));
    }

    // ======== associatedgraph: multiple polynomials ========
    CLPOLY_TEST("associatedgraph_multi");
    {
        variable x1("am1"), x2("am2"), x3("am3");
        polynomial_ZZ f1 = x1*x2;     // connects x1-x2
        polynomial_ZZ f2 = x2*x3;     // connects x2-x3
        std::vector<polynomial_ZZ> F = {f1, f2};
        auto G = associatedgraph(F);

        CLPOLY_ASSERT_EQ((int)G.size(), 3);
        CLPOLY_ASSERT_TRUE(G.is_edge(x1, x2));
        CLPOLY_ASSERT_TRUE(G.is_edge(x2, x3));
        CLPOLY_ASSERT_FALSE(G.is_edge(x1, x3));
    }

    // ======== peo: perfect elimination ordering ========
    CLPOLY_TEST("peo_basic");
    {
        variable x1("pe1"), x2("pe2"), x3("pe3"), x4("pe4");
        // Create a chain: x[i]*x[i+3] - x[i+1]*x[i+2]
        std::vector<polynomial_ZZ> F;
        F.push_back(x1*x4 - x2*x3);
        auto p = peo(F);
        // peo should return a permutation of all variables
        CLPOLY_ASSERT_EQ((int)p.size(), 4);
    }

    // ======== chordal_completion ========
    CLPOLY_TEST("chordal_completion_basic");
    {
        variable x1("cc1"), x2("cc2"), x3("cc3"), x4("cc4");
        polynomial_ZZ f1 = x1*x4 - x2*x3;
        std::vector<polynomial_ZZ> F = {f1};
        auto G = associatedgraph(F);
        auto Gc = chordal_completion(G);

        // Chordal completion should have at least as many edges
        CLPOLY_ASSERT(Gc.size() >= G.size());
    }

    // ======== connected_branch_graph ========
    CLPOLY_TEST("connected_branch_graph_basic");
    {
        variable x1("cb1"), x2("cb2"), x3("cb3");
        polynomial_ZZ f1 = x1*x2;
        polynomial_ZZ f2 = x2*x3;
        std::vector<polynomial_ZZ> F = {f1, f2};
        auto G = associatedgraph(F);
        auto Gcb = connected_branch_graph(G);

        // After connecting branches, x1-x3 should also be connected
        CLPOLY_ASSERT_TRUE(Gcb.is_edge(x1, x3));
    }

    // ======== graph_diff_score ========
    CLPOLY_TEST("graph_diff_score_basic");
    {
        variable x1("ds1"), x2("ds2"), x3("ds3");
        polynomial_ZZ f1 = x1*x2;
        polynomial_ZZ f2 = x2*x3;
        std::vector<polynomial_ZZ> F = {f1, f2};
        auto G = associatedgraph(F);

        // Same graph should have diff_score 0
        auto score_same = graph_diff_score(G, G);
        CLPOLY_ASSERT_EQ(score_same, 0.0);

        // Chordal completion should have diff_score >= 0
        auto Gc = chordal_completion(G);
        auto score_diff = graph_diff_score(Gc, G);
        CLPOLY_ASSERT(score_diff >= 0.0);
    }

    // ======== elimination_game ========
    CLPOLY_TEST("elimination_game_basic");
    {
        variable x1("eg1"), x2("eg2"), x3("eg3"), x4("eg4");
        std::vector<polynomial_ZZ> F;
        F.push_back(x1*x4 - x2*x3);
        auto G = associatedgraph(F);
        auto p = peo(F);
        auto Ge = elimination_game(G, p);
        // elimination_game should produce a graph with same number of nodes
        CLPOLY_ASSERT_EQ((int)Ge.size(), (int)G.size());
    }

    // ======== graph accessors ========
    CLPOLY_TEST("graph_accessors");
    {
        graph<variable> G;
        variable x1("ga1"), x2("ga2"), x3("ga3");
        G.add_node(x1);
        G.add_node(x2);
        G.add_node(x3);
        G.add_edge(x1, x2);

        // nodes()
        auto& nodes = G.nodes();
        CLPOLY_ASSERT_EQ((int)nodes.size(), 3);

        // adjacency_list()
        auto& adj = G.adjacency_list();
        CLPOLY_ASSERT_EQ((int)adj.size(), 3);

        // node(idx) and index(node)
        for (size_t i = 0; i < nodes.size(); ++i) {
            CLPOLY_ASSERT_EQ(G.node(i), nodes[i]);
            CLPOLY_ASSERT_EQ(G.index(nodes[i]), (uint64_t)i);
        }
    }

    // ======== graph index-based operations ========
    CLPOLY_TEST("graph_index_operations");
    {
        graph<variable> G;
        variable x1("gi1"), x2("gi2"), x3("gi3");
        G.add_node(x1);
        G.add_node(x2);
        G.add_node(x3);

        auto i1 = G.index(x1);
        auto i2 = G.index(x2);
        auto i3 = G.index(x3);

        // add_edge by index
        G.add_edge(i1, i2);
        CLPOLY_ASSERT_TRUE(G.is_edge(i1, i2));
        CLPOLY_ASSERT_TRUE(G.is_edge(i2, i1));  // undirected
        CLPOLY_ASSERT_FALSE(G.is_edge(i1, i3));

        // Also verify by node
        CLPOLY_ASSERT_TRUE(G.is_edge(x1, x2));
        CLPOLY_ASSERT_FALSE(G.is_edge(x1, x3));

        // Add another edge
        G.add_edge(i2, i3);
        CLPOLY_ASSERT_TRUE(G.is_edge(i2, i3));
        CLPOLY_ASSERT_TRUE(G.is_edge(x2, x3));
    }

    // ======== graph output ========
    CLPOLY_TEST("graph_output");
    {
        graph<variable> G;
        variable x1("go1"), x2("go2");
        G.add_node(x1);
        G.add_node(x2);
        G.add_edge(x1, x2);

        std::ostringstream oss;
        oss << G;
        std::string output = oss.str();
        CLPOLY_ASSERT_FALSE(output.empty());
    }

    // ======== elimination_height ========
    CLPOLY_TEST("graph_elimination_height");
    {
        variable x1("eh1"), x2("eh2"), x3("eh3"), x4("eh4");
        polynomial_ZZ f1 = x1*x2 + x3;
        polynomial_ZZ f2 = x2*x3 + x4;
        std::vector<polynomial_ZZ> F = {f1, f2};
        auto G = associatedgraph(F);
        auto p = peo(F);
        auto h = elimination_height(G, p);
        // Height should be at least 1 for a non-trivial graph
        CLPOLY_ASSERT_TRUE(h >= 1);
    }

    // ======== connected_branch ========
    CLPOLY_TEST("graph_connected_branch");
    {
        graph<variable> G;
        variable x1("cb_1"), x2("cb_2"), x3("cb_3"), x4("cb_4");
        G.add_node(x1);
        G.add_node(x2);
        G.add_node(x3);
        G.add_node(x4);
        // Two connected components: {x1, x2} and {x3, x4}
        G.add_edge(x1, x2);
        G.add_edge(x3, x4);

        auto branches = connected_branch(G);
        CLPOLY_ASSERT_EQ((int)branches.size(), 2);
        // Each branch should have 2 nodes
        for (auto& branch : branches) {
            CLPOLY_ASSERT_EQ((int)branch.size(), 2);
        }

        // Fully connected graph: one branch
        G.add_edge(x2, x3);
        auto branches2 = connected_branch(G);
        CLPOLY_ASSERT_EQ((int)branches2.size(), 1);
        CLPOLY_ASSERT_EQ((int)branches2[0].size(), 4);
    }

    // ======== perfect_elimination_ordering ========
    CLPOLY_TEST("graph_perfect_elim_ordering");
    {
        graph<variable> G;
        variable x1("peo1"), x2("peo2"), x3("peo3");
        G.add_node(x1);
        G.add_node(x2);
        G.add_node(x3);
        // Complete graph: all pairs connected
        G.add_edge(x1, x2);
        G.add_edge(x1, x3);
        G.add_edge(x2, x3);

        auto ordering = perfect_elimination_ordering(G);
        CLPOLY_ASSERT_EQ((int)ordering.size(), 3);
        // Should contain all nodes
        std::set<variable> nodes_set(ordering.begin(), ordering.end());
        CLPOLY_ASSERT_EQ((int)nodes_set.size(), 3);
        CLPOLY_ASSERT_TRUE(nodes_set.count(x1) > 0);
        CLPOLY_ASSERT_TRUE(nodes_set.count(x2) > 0);
        CLPOLY_ASSERT_TRUE(nodes_set.count(x3) > 0);
    }

    return clpoly_test::test_summary();
}
