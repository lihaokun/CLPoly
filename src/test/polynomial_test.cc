#include "polynomial.hh"
#include "variable.hh"
#include "monomial.hh"
#include <iostream>
int main(){
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::polynomial_ZZ p=2*x+y+1;
    std::cout<< p<<std::endl;
    auto l=clpoly::get_var(p);
    for (auto &i:l)
        std::cout<<i.first<<":"<<i.second<<" ";
    std::cout<<std::endl<<clpoly::get_deg(p)<<std::endl;
    
    return 0;
}