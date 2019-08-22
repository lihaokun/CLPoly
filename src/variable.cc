/*
Module Name:
    variable.cc
Abstract:
    定义类：variable
Author:
    haokun li
Notes:
*/
#include <iostream>
#include <unordered_map>
#include <vector>
#include "variable.hh"

namespace clpoly{
    std::vector<std::string> variable::variables={""};
    std::unordered_map<std::string,std::size_t> variable::name_map={{"",0}};
    std::size_t variable::init_variable(const std::string & variable_name)
    {
        #ifdef DEBUG
            if (variable::name_map.find(variable_name)!=variable::name_map.end())
                std::cout<<"Warning: Duplicate variable name.";
        #endif
        variable::name_map[variable_name]=variable::variables.size();
        variable::variables.push_back(variable_name);
        return variable::variables.size()-1;
        
    }
    variable variable::get_variable(const std::string & variable_name)
    {
        auto tmp=variable::name_map.find(variable_name);
        if (tmp==variable::name_map.end())
        {
            return variable(variable::init_variable(variable_name));
        }
        return variable(tmp->second);
    }
    variable variable::get_variable(std::size_t serial)
    {
        if (serial>=variable::variables.size()){
            throw std::invalid_argument("Serial of var not exist.");
        }
        return variable(serial);
    }
    variable variable::new_variable()
    {
        return variable(variable::init_variable("var_"+ std::to_string(variable::variables.size())));
    }
    variable variable::new_variable(const std::string & variable_name)
    {
        return variable(variable::init_variable(variable_name));
    }
}