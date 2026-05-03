// B2B C++ 端函数注册表
//
// dispatch(fn_name, args_json_array) → result_json
// 找不到函数抛 std::runtime_error("unknown fn: <name>")
//
// 添加新函数：在 cpp_registry.cc 加一个 if-分支

#ifndef CLPOLY_B2B_REGISTRY_HH
#define CLPOLY_B2B_REGISTRY_HH

#include <string>
#include "nlohmann/json.hpp"

namespace b2b {

nlohmann::json dispatch(const std::string& fn, const nlohmann::json& args);

}  // namespace b2b

#endif
