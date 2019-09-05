/*
Module Name:
    monomial.cc
Abstract:
    定义类：monomial
Author:
    haokun li
Notes:
*/
#include "monomial.hh"
#include "variable.hh"
#include <functional>
#include <sstream>
namespace clpoly{
    const std::function<bool(const variable &,const variable &)> monomial::init_comp(less<variable>);
    std::string  monomial::str() const {
        if (this->empty())
            return "1";
        std::ostringstream ss;
        bool is_print=false;
        for (auto&&i:this->__data)
        {
            if (i.second!=0)
            {
                if (is_print)
                    ss<<"*";
                ss<<i.first;
                is_print=true;
                if (i.second!=1)
                    ss<<"^"<<i.second;
            }
        }
        return ss.str();
    }
}