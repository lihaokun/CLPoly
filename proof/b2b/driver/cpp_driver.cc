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

using nlohmann::json;

int main() {
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
