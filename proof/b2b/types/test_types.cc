// 类型库 round-trip 单元测试：parse(serialize(x)) == x（按值语义比较）

#include "b2b_types.hh"
#include <iostream>
#include <cassert>
#include <vector>

using namespace b2b;

#define CHECK(cond, label) do { \
    if (!(cond)) { std::cerr << "FAIL: " << label << "\n"; return 1; } \
    else { std::cout << "  PASS: " << label << "\n"; } \
} while (0)

int main() {
    std::cout << "=== B2B types round-trip ===\n";

    // Int32
    {
        std::vector<int32_t> vs = {0, 1, -1, 2147483647, -2147483648};
        for (int32_t v : vs) {
            auto j = serialize_Int32(v);
            CHECK(parse_Int32(j) == v, "Int32 " + std::to_string(v));
        }
    }
    // Int64
    {
        std::vector<int64_t> vs = {0, 1, -1,
                                    9223372036854775807LL,
                                    -9223372036854775807LL - 1};
        for (int64_t v : vs) {
            auto j = serialize_Int64(v);
            CHECK(parse_Int64(j) == v, "Int64 " + std::to_string(v));
        }
    }
    // UInt64
    {
        std::vector<uint64_t> vs = {0u, 1u, 13u, 18446744073709551615ULL};
        for (uint64_t v : vs) {
            auto j = serialize_UInt64(v);
            CHECK(parse_UInt64(j) == v, "UInt64 " + std::to_string(v));
        }
    }
    // ZZ — small + large
    {
        std::vector<ZZ> zs;
        zs.emplace_back(0);
        zs.emplace_back(1);
        zs.emplace_back(-1);
        zs.emplace_back(123456789012345LL);
        zs.emplace_back("99999999999999999999999999");  // 26 digits
        for (const auto& z : zs) {
            auto j = serialize_ZZ(z);
            CHECK(parse_ZZ(j) == z, std::string("ZZ ") + j["val"].get<std::string>());
        }
    }
    // Zp
    {
        std::vector<Zp> zs = {Zp((int64_t)3, (uint64_t)7),
                              Zp((int64_t)0, (uint64_t)13),
                              Zp((int64_t)6, (uint64_t)7)};
        for (const auto& z : zs) {
            auto j = serialize_Zp(z);
            Zp parsed = parse_Zp(j);
            CHECK(parsed.number() == z.number() && parsed.prime() == z.prime(),
                  "Zp(" + std::to_string(z.number()) + "," + std::to_string(z.prime()) + ")");
        }
    }
    // SparsePolyZp
    {
        upolynomial_<Zp> p1;  // empty
        upolynomial_<Zp> p2;
        p2.push_back({umonomial(2), Zp((int64_t)3, (uint64_t)7)});
        p2.push_back({umonomial(0), Zp((int64_t)5, (uint64_t)7)});
        std::vector<upolynomial_<Zp>> ps = {p1, p2};
        for (const auto& p : ps) {
            auto j = serialize_SparsePolyZp(p);
            auto parsed = parse_SparsePolyZp(j);
            bool ok = parsed.size() == p.size();
            for (size_t i = 0; ok && i < p.size(); ++i) {
                ok = parsed[i].first.deg() == p[i].first.deg()
                  && parsed[i].second.number() == p[i].second.number()
                  && parsed[i].second.prime() == p[i].second.prime();
            }
            CHECK(ok, "SparsePolyZp size=" + std::to_string(p.size()));
        }
    }

    // 类型不匹配检测
    {
        bool threw = false;
        try { parse_Int32(serialize_Int64(7)); } catch (const std::exception&) { threw = true; }
        CHECK(threw, "type mismatch raises");
    }

    std::cout << "\n=== All round-trip tests PASSED ===\n";
    return 0;
}
