// B2B C++ driver — stdin/stdout NDJSON dispatcher
//
// 协议：
//   请求：{"id":"<str>", "fn":"<name>", "args":[...]}
//   响应：{"id":"<str>", "ok":true,  "ret":<value>}
//      或 {"id":"<str>", "ok":false, "err":"<msg>"}

#include <iostream>
#include <string>
#include "nlohmann/json.hpp"
#include "../registry/cpp_registry.hh"
#include "clpoly/variable.hh"

using nlohmann::json;

int main() {
    // 预声明 B2B 测试用变量，让 get_variable(serial) 调用合法。
    // Convention: serial 1=x, 2=y, 3=z, 4=w（与 vectors 中 Variable val 对齐）
    clpoly::variable _x("x"), _y("y"), _z("z"), _w("w");
    (void)_x; (void)_y; (void)_z; (void)_w;

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        json resp;
        std::string id;
        try {
            json req = json::parse(line);
            id = req.value("id", "");
            resp["id"] = id;
            std::string fn = req.at("fn").get<std::string>();
            json ret = b2b::dispatch(fn, req.at("args"));
            resp["ok"] = true;
            resp["ret"] = ret;
        } catch (const std::exception& e) {
            resp["id"] = id;
            resp["ok"] = false;
            resp["err"] = e.what();
        }
        std::cout << resp.dump() << "\n" << std::flush;
    }
    return 0;
}
