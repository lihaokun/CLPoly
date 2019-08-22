/*
Module Name:
    atomic_polynomial.hh
Abstract:
    定义类：atomic_polynomial
Author:
    haokun li
Notes:
*/
#include "variable.hh"
#include <vector>
namespace clpoly{
    template <typename Tm,typename Tc>
    class atomic_polynomial: public std::vector<std::pair<Tc,Tm>>
    {
        public:
            using std::vector<std::pair<Tc,Tm>>::vector;

    };
    

}
