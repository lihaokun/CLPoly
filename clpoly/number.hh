/*
Module Name:
    number.hh
Abstract:
    关于高精度的定义
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_NUMBER_HH
#define CLPOLY_NUMBER_HH
#include <gmpxx.h>
#include "basic.hh"
// #include <cln/integer.h>
// #include <cln/integer_io.h>
namespace clpoly{
    typedef mpz_class ZZ;
    typedef mpq_class QQ;
    
    // using ZZ=cln::cl_I;
    template<>
    inline void addmul(mpz_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_addmul(op.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void submul(mpz_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_submul(op.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void __div(mpz_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_fdiv_q(op.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void __div(mpq_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        op=mpq_class(op1,op2);
    }
    template <>
    void set_zero(mpz_class& op)
    {
        mpz_set_si(op.get_mpz_t(),0);
    }
    template<>
    struct zore_check<mpz_class>: public std::unary_function<mpz_class, bool>
    {
        bool operator()(const mpz_class & op)
        {
            return !op;
        } 
    };
}
#endif