#include "clpoly.hh"

int main(int argc, char const *argv[])
{
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::polynomial_ZZ f=-6*pow(x,4)-7*pow(x,2)*y-x*pow(y,3)*z-3*x*y-9;
    clpoly::polynomial_ZZ g=10*pow(x,4)*y+2*pow(x,4)-7*pow(x,3)*pow(y,2)+9*pow(x,2)+9*pow(z,4)-9*y+2*z;
    std::cout<<"f="<<f<<std::endl;
    std::cout<<"g="<<g<<std::endl;
    std::cout<<"prem(f,g,x)="<<prem(f,g,x)<<std::endl;
}
