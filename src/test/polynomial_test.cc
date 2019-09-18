#include "polynomial.hh"
#include "variable.hh"
#include "monomial.hh"
#include <iostream>
int main(){
    clpoly::variable x("x");
    clpoly::monomial m=x;
    clpoly::polynomial p=std::move(m);
    std::cout<<m<<std::endl;
    std::cout<<p<<std::endl;
    return 0;
}