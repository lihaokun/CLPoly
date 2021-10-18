/**
 * @file variable.cc
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义类：variable
 * 
 * 
 */
#include <clpoly/variable.hh>
namespace clpoly{
    std::vector<std::string> variable::variables={""};
    std::unordered_map<std::string,std::size_t> variable::name_map={{"",0}};
    std::vector<std::size_t>  variable::free_serial;
}