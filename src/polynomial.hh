/*
Module Name:
    atomic_polynomial.hh
Abstract:
    定义类：polynomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_POLYNOMIAL_HH
#define CLPOLY_POLYNOMIAL_HH
#include "basic.hh"
#include "variable.hh"
#include "monomial.hh"
#include "number.hh"
#include "atomic_polynomial.hh"
#include <vector>
#include <string>
namespace clpoly{
    class polynomial
    {
        private:
            atomic_polynomial<monomial,integer> __data;
        public:
            polynomial(){}
            polynomial(int64_t c):__data({{monomial(),c}}){}
            polynomial(const variable& v):__data({{monomial(v),1}}){}
            polynomial(const monomial & m):__data({{m,1}}){}
            polynomial(monomial && m):__data({{std::move(m),1}}){}
            polynomial(atomic_polynomial<monomial,integer> && p):__data(p){}
            polynomial(const polynomial & p):__data(p.__data){}
            polynomial(polynomial && p):__data(std::move(p.__data)){}
            std::string str() const
            {
                return this->__data.str();
            }
            friend std::ostream& operator<<  (std::ostream& stream, const polynomial & p)
            {
                return stream<<(p.str());
            }  
            
    };

}
#endif