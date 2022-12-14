#include <iostream>
#include <vector>
#include <clpoly/clpoly.hh>
#include <functional>

int main(int argc, char const *argv[])
{
    clpoly::interval I1;
    // std::cout<<I1;
    I1.set_l(1);
    I1.set_r(10);
    std::cout<<"I1:="<<I1<<std::endl;
    clpoly::interval I2;
    // std::cout<<I2;
    I2.set_l(2);
    I2.set_r(20);
    std::cout<<"I2:="<<I2<<std::endl;
    std::cout<<"I1+I2=="<<I1+I2<<std::endl;
    std::cout<<"I1-I2=="<<I1-I2<<std::endl;
    std::cout<<"I1*I2=="<<I1*I2<<std::endl;
    std::cout<<"I1/I2=="<<I1/I2<<std::endl;
    std::cout<<"(I1-I2)^3=="<<(I1-I2).pow(3)<<std::endl;

    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    
    clpoly::polynomial_ZZ  p=y*x+pow(y,3);
    std::cout<<"p="<<p<<std::endl;
    std::cout<<"p="<<clpoly::assign(p,{std::make_pair(y,I1)})<<std::endl;
    
    return 0;
}
