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
#include <functional>
namespace clpoly{
    const std::function<bool(const variable &,const variable &)> monomial::init_comp(less<variable>);
}