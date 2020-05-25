/*
Module Name:
    monomial.hh
Abstract:
    定义类：monomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_MONOMIAL_HH
#define CLPOLY_MONOMIAL_HH
#include "variable.hh"
#include "basic.hh"
#include <vector>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include "monomial_order.hh"
#include "basic_monomial.hh"

namespace clpoly
{
    using monomial=basic_monomial<grlex>;
}

#endif