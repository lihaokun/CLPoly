// B2B C++ 函数注册表实现
//
// 每加一个函数：写一个 if (fn == "<name>") { ... return ...; } 分支

#include "cpp_registry.hh"
#include "../types/b2b_types.hh"

#include "clpoly/polynomial_factorize_zp.hh"
#include "clpoly/polynomial_factorize_univar.hh"
#include "clpoly/polynomial.hh"
#include <cstring>
#include <gmp.h>
#include <sstream>

namespace b2b {

using nlohmann::json;

json dispatch(const std::string& fn, const json& args) {
    if (!args.is_array()) {
        throw std::runtime_error("args: expected array");
    }

    // === 函数 dispatch ===
    if (fn == "__lll_factorize") {
        auto f = parse_SparsePolyZZ(args.at(0));
        std::vector<clpoly::upolynomial_<clpoly::Zp>> modular;
        for (const auto& poly : args.at(1)) modular.push_back(parse_SparsePolyZp(poly));
        auto prime = args.at(2).get<uint64_t>();
        if (prime < 2 || modular.size() < 2 || modular.size() >= 2147483648ULL || clpoly::get_deg(f) <= 0)
            throw std::runtime_error("invalid retry test input");
        auto result = clpoly::__lll_factorize(f, modular, prime);
        std::vector<json> factors;
        for (const auto& factor : result) factors.push_back(serialize_SparsePolyZZ(factor));
        std::sort(factors.begin(), factors.end(),
            [](const json& a, const json& b) { return a.dump() < b.dump(); });
        return factors;
    }
    if (fn == "__vanhoeij_try_recombine") {
        auto f = parse_SparsePolyZZ(args.at(0));
        std::vector<clpoly::upolynomial_<clpoly::ZZ>> lifted;
        for (const auto& poly : args.at(1)) lifted.push_back(parse_SparsePolyZZ(poly));
        auto modulus = parse_ZZ(args.at(2));
        if (lifted.size() < 2 || lifted.size() >= 2147483648ULL || modulus <= 1 || clpoly::get_deg(f) <= 0)
            throw std::runtime_error("invalid try-recombine test input");
        auto result = clpoly::__vanhoeij_try_recombine(f, lifted, modulus);
        std::vector<json> factors;
        for (const auto& factor : result.factors) factors.push_back(serialize_SparsePolyZZ(factor));
        auto less = [](const json& a, const json& b) { return a.dump() < b.dump(); };
        std::sort(factors.begin(), factors.end(), less);
        auto expected = args.at(4).get<std::vector<json>>();
        std::sort(expected.begin(), expected.end(), less);
        if (result.complete != args.at(3).get<bool>() || factors != expected)
            throw std::runtime_error("try-recombine differs from independently known factors/completion");
        return json::array({result.complete, factors, result.cld_columns, result.retained_rows});
    }
    if (fn == "__vanhoeij_choose_parameters") {
        int count = args.at(0).get<int>();
        auto height = parse_ZZ(args.at(1));
        auto modulus = parse_ZZ(args.at(2));
        if (count < 2 || height < 0 || modulus <= 1)
            throw std::runtime_error("require count>=2, height>=0, modulus>1");
        clpoly::__vanhoeij_parameters output{clpoly::ZZ(17), clpoly::ZZ(19), clpoly::ZZ(23)};
        if (!clpoly::__vanhoeij_choose_parameters(count, height, modulus, output)) {
            if (output.divisor != 17 || output.modulus != 19 || output.residual_bound != 23)
                throw std::runtime_error("failed parameter selection changed output");
            return json::array();
        }
        return json::array({serialize_ZZ(output.divisor), serialize_ZZ(output.modulus),
                            serialize_ZZ(output.residual_bound)});
    }
    if (fn == "__vanhoeij_cld_width") {
        std::vector<clpoly::upolynomial_<clpoly::ZZ>> cld;
        for (const auto& input : args.at(0)) {
            auto poly = parse_SparsePolyZZ(input);
            if (!poly.empty() && clpoly::get_deg(poly) >= 2147483647)
                throw std::runtime_error("degree+1 must fit Int32");
            cld.push_back(std::move(poly));
        }
        return clpoly::__vanhoeij_cld_width(cld);
    }
    if (fn == "__vanhoeij_singleton_groups") {
        int count = args.at(0).get<int>();
        if (count < 0 || count > 4096) throw std::runtime_error("test count must be in 0..4096");
        return clpoly::__vanhoeij_singleton_groups(count);
    }
    if (fn == "__vanhoeij_next_target") {
        int target = args.at(0).get<int>();
        int width = args.at(1).get<int>();
        if (target <= 0 || target > width) throw std::runtime_error("require 0 < target <= width");
        return clpoly::__vanhoeij_next_target(target, width);
    }
    if (fn == "__vanhoeij_clamp_precision") {
        const uint64_t bits = args.at(0).get<uint64_t>();
        double value;
        std::memcpy(&value, &bits, sizeof(value));
        return clpoly::__vanhoeij_clamp_precision(value);
    }
    if (fn == "__vanhoeij_full_precision") {
        const uint64_t prime = args.at(0).get<uint64_t>();
        if (prime < 2) throw std::runtime_error("require p >= 2");
        return clpoly::__vanhoeij_full_precision(prime, parse_ZZ(args.at(1)));
    }
    if (fn == "__heuristic_starting_precision") {
        auto f = parse_SparsePolyZZ(args.at(0));
        int count = args.at(1).get<int>();
        uint64_t prime = args.at(2).get<uint64_t>();
        if (f.empty() || count < 0 || prime < 2) throw std::runtime_error("invalid heuristic input");
        auto result = clpoly::__heuristic_starting_precision(f, count, prime);
        return json::array({result.first, result.second});
    }
    if (fn == "__vanhoeij_next_precision") {
        int precision = args.at(0).get<int>();
        int full = args.at(1).get<int>();
        if (precision < 0 || full < 0) throw std::runtime_error("require nonnegative precision");
        try { return json{{"value", clpoly::__vanhoeij_next_precision(precision, full)}}; }
        catch (const std::overflow_error&) { return json{{"overflow", true}}; }
    }
    if (fn == "__vanhoeij_validate_partition") {
        auto f = parse_SparsePolyZZ(args.at(0));
        std::vector<clpoly::upolynomial_<clpoly::ZZ>> lifted;
        for (const auto& poly : args.at(1)) lifted.push_back(parse_SparsePolyZZ(poly));
        auto modulus = parse_ZZ(args.at(2));
        if (modulus <= 1) throw std::runtime_error("modulus must exceed one");
        auto groups = args.at(3).get<std::vector<std::vector<int>>>();
        for (const auto& group : groups) for (int index : group)
            if (index < 0 || (size_t)index >= lifted.size())
                throw std::runtime_error("invalid physical group index");
        clpoly::__vanhoeij_result answer;
        clpoly::upolynomial_<clpoly::ZZ> sentinel;
        sentinel.push_back({clpoly::umonomial(0), clpoly::ZZ(17)});
        answer.factors = {sentinel};
        answer.complete = args.at(4).get<bool>();
        answer.cld_columns = 7;
        answer.retained_rows = 9;
        bool accepted = clpoly::__vanhoeij_validate_partition(f, lifted, modulus, groups, answer);
        std::vector<json> factors;
        for (const auto& poly : answer.factors) factors.push_back(serialize_SparsePolyZZ(poly));
        // Equal-degree std::sort ties are unspecified: compare the resulting
        // multiset while testing all state and acceptance fields exactly.
        std::sort(factors.begin(), factors.end(), [](const json& a, const json& b) {
            return a.dump() < b.dump();
        });
        return json::array({accepted, answer.complete, factors, answer.cld_columns, answer.retained_rows});
    }
    if (fn == "__vanhoeij_cld_bound") {
        auto f = parse_SparsePolyZZ(args.at(0));
        if (clpoly::get_deg(f) <= 0 || clpoly::get_deg(f) > 4096)
            throw std::runtime_error("test degree must be in 1..4096");
        return serialize_ZZ(clpoly::__vanhoeij_cld_bound(f));
    }
    if (fn == "__upoly_norm_l1") {
        return serialize_ZZ(clpoly::__upoly_norm_l1(parse_SparsePolyZZ(args.at(0))));
    }
    if (fn == "__vanhoeij_truncate_cld") {
        std::vector<clpoly::upolynomial_<clpoly::ZZ>> cld;
        for (const auto& input : args.at(0)) cld.push_back(parse_SparsePolyZZ(input));
        auto divisor = parse_ZZ(args.at(1));
        if (divisor <= 0) throw std::runtime_error("positive divisor required");
        clpoly::__vanhoeij_truncate_cld(cld, divisor);
        json output = json::array();
        for (const auto& poly : cld) output.push_back(serialize_SparsePolyZZ(poly));
        return output;
    }
    if (fn == "__lll_reduce_with_norms") {
        clpoly::LLLMatrix matrix;
        for (const auto& input_row : args.at(0)) {
            std::vector<clpoly::ZZ> row;
            for (const auto& entry : input_row) row.push_back(parse_ZZ(entry));
            matrix.push_back(std::move(row));
        }
        clpoly::LLLMatrix transform;
        std::vector<clpoly::QQ> norms;
        clpoly::__lll_reduce(matrix, transform, parse_ZZ(args.at(1)), &norms);
        json encoded_norms = json::array();
        for (const auto& norm : norms) {
            std::ostringstream numerator, denominator;
            numerator << norm.get_num();
            denominator << norm.get_den();
            encoded_norms.push_back(json::array({numerator.str(), denominator.str()}));
        }
        // The repaired caller observes the basis, U and GS cache. It ignores
        // the legacy individual-short-row return (whose equal-norm sort ties
        // are not specified by std::sort).
        auto encode_matrix = [](const clpoly::LLLMatrix& input) {
            json output = json::array();
            for (const auto& row : input) {
                json entries = json::array();
                for (const auto& entry : row) entries.push_back(serialize_ZZ(entry));
                output.push_back(entries);
            }
            return output;
        };
        return json::array({encode_matrix(matrix), encode_matrix(transform), encoded_norms});
    }
    if (fn == "__lll_short_span_from_norms") {
        std::vector<clpoly::QQ> norms;
        for (const auto& pair : args.at(0)) {
            auto denominator = parse_ZZ(pair.at(1));
            if (denominator == 0) throw std::runtime_error("zero denominator");
            norms.emplace_back(parse_ZZ(pair.at(0)), denominator);
        }
        json result = json::array();
        for (size_t row : clpoly::__lll_short_span_from_norms(norms, parse_ZZ(args.at(1))))
            result.push_back(serialize_UInt64(static_cast<uint64_t>(row)));
        return result;
    }
    if (fn == "__build_cld_matrix") {
        clpoly::LLLMatrix matrix;
        for (const auto& input_row : args.at(0)) {
            std::vector<clpoly::ZZ> row;
            for (const auto& entry : input_row) row.push_back(parse_ZZ(entry));
            matrix.push_back(std::move(row));
        }
        std::vector<clpoly::upolynomial_<clpoly::ZZ>> cld;
        for (const auto& input : args.at(1)) cld.push_back(parse_SparsePolyZZ(input));
        auto current = parse_UInt64(args.at(2));
        auto target = parse_UInt64(args.at(3));
        if (current > INT_MAX || target > INT_MAX ||
            cld.size() > INT_MAX - current || matrix.size() != cld.size() + current)
            throw std::runtime_error("invalid lattice dimensions");
        for (const auto& row : matrix)
            if (row.size() != matrix.size()) throw std::runtime_error("non-square lattice");
        int added = clpoly::__build_cld_matrix(matrix, cld, static_cast<int>(current),
            static_cast<int>(target), parse_ZZ(args.at(4)));
        json output = json::array();
        for (const auto& row : matrix) {
            json entries = json::array();
            for (const auto& entry : row) entries.push_back(serialize_ZZ(entry));
            output.push_back(entries);
        }
        return json::array({output, serialize_UInt64(static_cast<uint64_t>(added))});
    }
    if (fn == "__factor_Zp") {
        auto f = parse_SparsePolyZp(args.at(0));
        return serialize_FactorZpResult(clpoly::__factor_Zp(std::move(f)));
    }
    if (fn == "__factor_squarefree_primitive_ZZ") {
        auto f = parse_SparsePolyZZ(args.at(0));
        bool use_large_prime = parse_Bool(args.at(1));
        bool previous = clpoly::__g_use_large_prime;
        clpoly::__g_use_large_prime = use_large_prime;
        try {
            auto result = clpoly::__factor_squarefree_primitive_ZZ(f);
            clpoly::__g_use_large_prime = previous;
            if (args.size() > 2) {
                std::vector<json> factors;
                for (const auto& factor : result) factors.push_back(serialize_SparsePolyZZ(factor));
                auto expected = args.at(2).get<std::vector<json>>();
                auto less = [](const json& a, const json& b) { return a.dump() < b.dump(); };
                std::sort(factors.begin(), factors.end(), less);
                std::sort(expected.begin(), expected.end(), less);
                if (factors != expected)
                    throw std::runtime_error("outer FactorZZ differs from independently known factors");
            }
            return serialize_ArraySPZZ(result);
        } catch (...) {
            clpoly::__g_use_large_prime = previous;
            throw;
        }
    }

    if (fn == "__hensel_factor_count_fits") {
        uint64_t factor_count = parse_UInt64(args.at(0));
        return serialize_Bool(clpoly::__hensel_factor_count_fits(
            static_cast<size_t>(factor_count)));
    }

    if (fn == "__make_zp") {
        int64_t  val = parse_Int64(args.at(0));
        uint64_t p   = parse_UInt64(args.at(1));
        return serialize_Zp(clpoly::__make_zp(val, p));
    }

    // ---- Zp 算术（B2B wrapper）：operator 重载在 C++ 没有命名函数，B2B 起 ____xx_zp 名 ----
    if (fn == "__add_zp") {
        Zp a = parse_Zp(args.at(0));
        Zp b = parse_Zp(args.at(1));
        return serialize_Zp(a + b);
    }
    if (fn == "__sub_zp") {
        Zp a = parse_Zp(args.at(0));
        Zp b = parse_Zp(args.at(1));
        return serialize_Zp(a - b);
    }
    if (fn == "__mul_zp") {
        Zp a = parse_Zp(args.at(0));
        Zp b = parse_Zp(args.at(1));
        return serialize_Zp(a * b);
    }
    if (fn == "__neg_zp") {
        Zp a = parse_Zp(args.at(0));
        return serialize_Zp(-a);
    }
    if (fn == "__inv_zp") {
        Zp a = parse_Zp(args.at(0));
        return serialize_Zp(a.inv());
    }
    if (fn == "__pow_zp") {
        Zp       a = parse_Zp(args.at(0));
        int64_t  i = parse_Int64(args.at(1));
        return serialize_Zp(clpoly::pow(a, i));
    }

    // ---- ZZ 算术（C++ 端用 gcd/lcm/static fdiv_q/r 入口）----
    if (fn == "__gcd_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        return serialize_ZZ(clpoly::gcd(a, b));
    }
    if (fn == "__lcm_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        return serialize_ZZ(clpoly::lcm(a, b));
    }
    if (fn == "__fdiv_q_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        ZZ q;
        ZZ::fdiv_q(q, a, b);
        return serialize_ZZ(q);
    }
    if (fn == "__fdiv_r_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        ZZ r;
        ZZ::fdiv_r(r, a, b);
        return serialize_ZZ(r);
    }
    if (fn == "__fdiv_ui_zz") {
        ZZ       a = parse_ZZ(args.at(0));
        uint64_t b = parse_UInt64(args.at(1));
        return serialize_UInt64(a.fdiv_ui(b));
    }
    if (fn == "__sizeinbase_zz") {
        ZZ      z    = parse_ZZ(args.at(0));
        int32_t base = parse_Int32(args.at(1));
        return serialize_UInt64(z.sizeinbase(base));
    }
    if (fn == "__invert_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ m = parse_ZZ(args.at(1));
        ZZ result;
        bool ok = ZZ::invert(result, a, m);
        return serialize_BoolZZ(ok, result);
    }

    // ---- SparsePolyZZ ops ----
    if (fn == "__cont_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        return serialize_ZZ(clpoly::cont(p));
    }
    if (fn == "__pp_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        ZZ c = clpoly::cont(p);
        // pp = p / cont(p)（精确除法）。空 p 或 cont=0：返原值（防 div by zero）
        if (c != ZZ(0)) {
            for (auto& term : p) term.second = term.second / c;
        }
        return serialize_SparsePolyZZ(p);
    }
    if (fn == "__derivative_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        return serialize_SparsePolyZZ(clpoly::derivative(p));
    }
    if (fn == "__normalization_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        p.normalization();
        return serialize_SparsePolyZZ(p);
    }
    if (fn == "__polynomial_mod_zz") {
        auto     p     = parse_SparsePolyZZ(args.at(0));
        uint64_t prime = parse_UInt64(args.at(1));
        return serialize_SparsePolyZp(clpoly::polynomial_mod(p, prime));
    }

    // ---- SparsePolyZp ops ----
    if (fn == "__add_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(f + g);
    }
    if (fn == "__sub_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(f - g);
    }
    if (fn == "__mul_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(f * g);
    }
    if (fn == "__divmod_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        upolynomial_<Zp> q, r;
        clpoly::__upoly_divmod(q, r, f, g);
        return serialize_PairSPZp(q, r);
    }
    if (fn == "__gcd_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(clpoly::polynomial_GCD(f, g));
    }
    if (fn == "__gcd_eea_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        upolynomial_<Zp> s, t;
        auto gcd = clpoly::polynomial_GCD(f, g, s, t);
        return serialize_TripleSPZp(gcd, s, t);
    }
    if (fn == "__derivative_zp_poly") {
        auto p = parse_SparsePolyZp(args.at(0));
        return serialize_SparsePolyZp(clpoly::derivative(p));
    }
    if (fn == "__normalization_zp_poly") {
        auto p = parse_SparsePolyZp(args.at(0));
        p.normalization();
        return serialize_SparsePolyZp(p);
    }

    // ---- A5: MvPoly ops ----
    if (fn == "__assign_mv_zz") {
        auto p = parse_MvPolyZZ(args.at(0));
        auto v = parse_Variable(args.at(1));
        auto c = parse_ZZ(args.at(2));
        return serialize_MvPolyZZ(clpoly::assign(p, v, c));
    }
    if (fn == "__assign2_mv_zz") {
        auto p = parse_MvPolyZZ(args.at(0));
        auto m = parse_VarMapZZ(args.at(1));
        return serialize_MvPolyZZ(clpoly::assign(p, m));
    }
    if (fn == "__leadcoeff_mv_zz") {
        auto p = parse_MvPolyZZ(args.at(0));
        auto v = parse_Variable(args.at(1));
        return serialize_MvPolyZZ(clpoly::leadcoeff(p, v));
    }
    if (fn == "__poly_convert_spzp_to_spzz") {
        auto p_in = parse_SparsePolyZp(args.at(0));
        upolynomial_<ZZ> p_out;
        clpoly::poly_convert(p_in, p_out);
        return serialize_SparsePolyZZ(p_out);
    }
    if (fn == "__poly_convert3_spzz_to_mvzz") {
        auto p_in = parse_SparsePolyZZ(args.at(0));
        auto v    = parse_Variable(args.at(1));
        MvPolyZZ p_out;
        clpoly::poly_convert(p_in, p_out, v);
        return serialize_MvPolyZZ(p_out);
    }

    // ---- A6: Bezout extGcd ----
    if (fn == "__nat_extgcd") {
        // 输入两个非负整数 a, b，返回 (g, s, t) 满足 s*a + t*b = g (= gcd(a,b))
        // 直接调 GMP mpz_gcdext，然后用字符串往返回 ZZ（绕过私有 mpz 访问）
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        mpz_t g_mp, s_mp, t_mp, a_mp, b_mp;
        mpz_inits(g_mp, s_mp, t_mp, NULL);
        // 把 ZZ 转 mpz_t（用 ZZ::operator<<）
        std::ostringstream a_os; a_os << a;
        std::ostringstream b_os; b_os << b;
        mpz_init_set_str(a_mp, a_os.str().c_str(), 10);
        mpz_init_set_str(b_mp, b_os.str().c_str(), 10);
        mpz_gcdext(g_mp, s_mp, t_mp, a_mp, b_mp);
        // 从 mpz_t 重建 ZZ
        char* g_str = mpz_get_str(nullptr, 10, g_mp);
        char* s_str = mpz_get_str(nullptr, 10, s_mp);
        char* t_str = mpz_get_str(nullptr, 10, t_mp);
        ZZ g(g_str), s(s_str), t(t_str);
        void (*freefunc)(void*, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(g_str, std::strlen(g_str) + 1);
        freefunc(s_str, std::strlen(s_str) + 1);
        freefunc(t_str, std::strlen(t_str) + 1);
        mpz_clears(g_mp, s_mp, t_mp, a_mp, b_mp, NULL);
        return serialize_ZZTriple(g, s, t);
    }

    throw std::runtime_error("unknown fn: " + fn);
}

}  // namespace b2b
