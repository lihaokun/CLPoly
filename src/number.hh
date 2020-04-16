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
namespace clpoly{
    typedef mpz_class ZZ;
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
}
#endif