#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <set>
#include <algorithm>

using namespace clpoly;

// Helper: count undirected edges in a graph via adjacency list
template<class T>
size_t count_edges(const graph<T>& G) {
    size_t total = 0;
    for (size_t i = 0; i < G.size(); ++i)
        total += G.adjacency_list()[i].size();
    return total / 2; // each undirected edge counted twice
}

int main() {
    // ================================================================
    //  random_select
    // ================================================================

    CLPOLY_TEST("random_select_basic");
    {
        auto v = random_select(100, 10);
        CLPOLY_ASSERT_EQ(v.size(), (size_t)10);
        for (auto x : v) {
            CLPOLY_ASSERT_TRUE(x < 100);
        }
        // Should be sorted
        for (size_t i = 1; i < v.size(); ++i) {
            CLPOLY_ASSERT_TRUE(v[i] >= v[i - 1]);
        }
    }

    CLPOLY_TEST("random_select_no_repeat");
    {
        auto v = random_select(20, 10, true);
        CLPOLY_ASSERT_EQ(v.size(), (size_t)10);
        // Check no duplicates
        std::set<size_t> s(v.begin(), v.end());
        CLPOLY_ASSERT_EQ(s.size(), v.size());
        // All values in range
        for (auto x : v) {
            CLPOLY_ASSERT_TRUE(x < 20);
        }
    }

    CLPOLY_TEST("random_select_n_ge_m");
    {
        auto v = random_select(5, 10, true);
        CLPOLY_ASSERT_EQ(v.size(), (size_t)5);
        for (size_t i = 0; i < 5; ++i) {
            CLPOLY_ASSERT_EQ(v[i], i);
        }
    }

    // ================================================================
    //  random_polynomial (fixed length version)
    // ================================================================

    CLPOLY_TEST("random_polynomial_fixed_len");
    {
        variable x("rp_x"), y("rp_y");
        std::vector<variable> vars = {x, y};
        auto p = random_polynomial<ZZ>(vars, 5, 4, {-10, 10});
        CLPOLY_ASSERT_FALSE(p.empty());
        CLPOLY_ASSERT_TRUE(p.size() <= 4);
        CLPOLY_ASSERT_EQ(p.degree(), (int64_t)5);
    }

    CLPOLY_TEST("random_polynomial_fixed_len_coeff_range");
    {
        variable x("rpc_x"), y("rpc_y");
        std::vector<variable> vars = {x, y};
        auto p = random_polynomial<ZZ>(vars, 3, 5, {1, 100});
        for (auto& term : p) {
            int c = (int)term.second.get_si();
            CLPOLY_ASSERT_TRUE(c >= 1 && c <= 100);
        }
    }

    CLPOLY_TEST("random_polynomial_fixed_len_add_num");
    {
        variable x("rpn_x"), y("rpn_y");
        std::vector<variable> vars = {x, y};
        auto p = random_polynomial<ZZ>(vars, 3, 5, {1, 100}, true);
        // Check that the polynomial has a constant term (last term, empty monomial)
        CLPOLY_ASSERT_FALSE(p.empty());
        bool has_const = false;
        for (auto& term : p) {
            if (term.first.empty()) {
                has_const = true;
                break;
            }
        }
        CLPOLY_ASSERT_TRUE(has_const);
    }

    // ================================================================
    //  random_polynomial (probability version)
    // ================================================================

    CLPOLY_TEST("random_polynomial_prob_p1");
    {
        variable x("rpp_x"), y("rpp_y");
        std::vector<variable> vars = {x, y};
        auto p = random_polynomial<ZZ>(vars, 3, 1.0, 1, 100);
        // p=1.0 means all monomials included (dense)
        CLPOLY_ASSERT_FALSE(p.empty());
    }

    CLPOLY_TEST("random_polynomial_prob_p0");
    {
        variable x("rp0_x"), y("rp0_y");
        std::vector<variable> vars = {x, y};
        auto p = random_polynomial<ZZ>(vars, 3, 0.0, 1, 100);
        // p=0.0 should give empty polynomial
        CLPOLY_ASSERT_TRUE(p.empty());
    }

    CLPOLY_TEST("random_polynomial_prob_degree");
    {
        variable x("rpd_x"), y("rpd_y");
        std::vector<variable> vars = {x, y};
        auto p = random_polynomial<ZZ>(vars, 4, 0.5, -10, 10);
        if (!p.empty()) {
            CLPOLY_ASSERT_TRUE(p.degree() <= 4);
        }
    }

    // ================================================================
    //  RandomSample
    // ================================================================

    CLPOLY_TEST("RandomSample_basic");
    {
        std::vector<int> input = {10, 20, 30, 40, 50};
        std::set<int> input_set(input.begin(), input.end());

        // n < size: returns n elements
        auto output = RandomSample(input, (uint64_t)3);
        CLPOLY_ASSERT_EQ(output.size(), (size_t)3);
        for (auto x : output) {
            CLPOLY_ASSERT_TRUE(input_set.count(x) > 0);
        }
        // No duplicates
        std::set<int> out_set(output.begin(), output.end());
        CLPOLY_ASSERT_EQ(out_set.size(), output.size());

        // n == size: returns all elements
        auto output2 = RandomSample(input, (uint64_t)5);
        CLPOLY_ASSERT_EQ(output2.size(), (size_t)5);

        // n > size: clamped to size
        auto output3 = RandomSample(input, (uint64_t)10);
        CLPOLY_ASSERT_EQ(output3.size(), (size_t)5);
    }

    // ================================================================
    //  random_graph
    // ================================================================

    CLPOLY_TEST("random_graph_p0");
    {
        std::vector<variable> nodes;
        for (int i = 0; i < 5; ++i)
            nodes.push_back(variable("rg0_" + std::to_string(i)));
        auto G = random_graph<variable>(nodes, 0.0);
        CLPOLY_ASSERT_EQ(G.size(), (size_t)5);
        CLPOLY_ASSERT_EQ(count_edges(G), (size_t)0);
    }

    CLPOLY_TEST("random_graph_p1");
    {
        std::vector<variable> nodes;
        int n = 5;
        for (int i = 0; i < n; ++i)
            nodes.push_back(variable("rg1_" + std::to_string(i)));
        auto G = random_graph<variable>(nodes, 1.0);
        CLPOLY_ASSERT_EQ(G.size(), (size_t)n);
        CLPOLY_ASSERT_EQ(count_edges(G), (size_t)(n * (n - 1) / 2));
    }

    // ================================================================
    //  random_polynomials
    // ================================================================

    CLPOLY_TEST("random_polynomials_basic");
    {
        std::vector<variable> vars;
        for (int i = 0; i < 4; ++i)
            vars.push_back(variable("rps_" + std::to_string(i)));
        auto polys = random_polynomials<ZZ>(vars, 3, 1.0, 0.8, 10, -10);
        // With p1=1.0 (complete graph), should produce some polynomials
        CLPOLY_ASSERT_FALSE(polys.empty());
        for (auto& p : polys) {
            if (!p.empty()) {
                CLPOLY_ASSERT_TRUE(p.degree() <= 3);
            }
        }
    }

    // ================================================================
    //  random_polynomial: single variable
    // ================================================================

    CLPOLY_TEST("random_polynomial_single_var");
    {
        variable x("rpsv_x");
        std::vector<variable> vars = {x};
        // With single variable, degree d => at most d+1 monomials
        auto p = random_polynomial<ZZ>(vars, 4, 3, {-10, 10});
        CLPOLY_ASSERT_FALSE(p.empty());
        CLPOLY_ASSERT_TRUE(p.degree() <= 4);
        CLPOLY_ASSERT_TRUE(p.size() <= 3);

        // Probability version with single variable
        auto p2 = random_polynomial<ZZ>(vars, 3, 1.0, 1, 100);
        CLPOLY_ASSERT_FALSE(p2.empty());
        CLPOLY_ASSERT_TRUE(p2.degree() <= 3);
    }

    // ================================================================
    //  random_polynomial: len clamping
    // ================================================================

    CLPOLY_TEST("random_polynomial_len_clamp");
    {
        variable x("rplc_x"), y("rplc_y");
        std::vector<variable> vars = {x, y};
        // degree=2, 2 vars => total monomials = C(2+2,2) = 6
        // Request len=100 (much more than possible)
        auto p = random_polynomial<ZZ>(vars, 2, 100, {1, 10});
        // Should clamp: at most 6 monomials
        CLPOLY_ASSERT_TRUE(p.size() <= 6);
        CLPOLY_ASSERT_TRUE(p.degree() <= 2);

        // degree=1, 2 vars => monomials: 1, x, y => 3
        auto p2 = random_polynomial<ZZ>(vars, 1, 50, {1, 10});
        CLPOLY_ASSERT_TRUE(p2.size() <= 3);
    }

    // ================================================================
    //  random_upolynomial
    // ================================================================

    CLPOLY_TEST("random_upolynomial_basic");
    {
        auto p = random_upolynomial<ZZ>(5, 4, {-10, 10});
        CLPOLY_ASSERT_EQ(degree(p), (int64_t)5);
        CLPOLY_ASSERT_TRUE(p.size() >= 1 && p.size() <= 6);
        for (auto& t : p)
            CLPOLY_ASSERT_TRUE(bool(t.second));
    }

    // ================================================================
    //  random_polynomial_QQ
    // ================================================================

    CLPOLY_TEST("random_polynomial_QQ_basic");
    {
        variable x("rpqq_x"), y("rpqq_y");
        auto p = random_polynomial_QQ({x, y}, 3, 4, {-10, 10}, 8);
        CLPOLY_ASSERT_EQ(degree(p), (int64_t)3);
        CLPOLY_ASSERT_TRUE(p.size() >= 1 && p.size() <= 5);
        for (auto& t : p) {
            CLPOLY_ASSERT_TRUE(bool(t.second));
            CLPOLY_ASSERT_TRUE(t.second.den() > ZZ(0));
        }
    }

    return clpoly_test::test_summary();
}
