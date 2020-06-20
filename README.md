# CLPoly

CLPply 是一个开发中的多项式c++库，目标是一个高效易且的多项式库。

CLPoly is a c++ library of polynomial.

## Example

```
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
```
结果
```
f=-x*y^3*z-6*x^4-7*x^2*y-3*x*y-9
g=10*x^4*y-7*x^3*y^2+2*x^4+9*z^4+9*x^2-9*y+2*z
prem(f,g,x)=-10*x*y^4*z-42*x^3*y^2-2*x*y^3*z-70*x^2*y^2+54*z^4-14*x^2*y-30*x*y^2+54*x^2-6*x*y-144*y+12*z-18
```
