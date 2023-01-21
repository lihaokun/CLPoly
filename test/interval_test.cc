#include <iostream>
#include <vector>
#include <clpoly/clpoly.hh>
#include <clpoly/interval.hh>
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
    // {
    //     clpoly::interval I(0,10);
    //     std::cout<<I<<std::endl;
    //     std::cout<<bool(I)<<std::endl;
    //     std::cout<<(I==0)<<std::endl;
    //     std::cout<<(I==1)<<std::endl;
    //     std::cout<<(I!=1)<<std::endl;
        
    // }

    clpoly::polynomial_ZZ  p=-y*pow(x,1)+2*y*pow(x,2)+pow(y,3)+3*z*y*pow(x,3);
    std::cout<<"p="<<p<<std::endl;
    {
        std::map<clpoly::variable,clpoly::interval> tmp;
        tmp[y]=I1;tmp[z]=I1;
        clpoly::interval_poly<> p1=clpoly::assign(p,tmp);
        std::cout<<"p="<<p1<<std::endl;
        // {
        // clpoly::QQ l1=p1.back().second.get_l();
        // clpoly::QQ l2=p1.front().second.get_r();
        // std::cout<<"l1/l2="<<l1<<" "<<l2<<" "<<l1/l2<<std::endl;
        // } 
        std::cout<< feasible_range(p1,'=') <<std::endl;
        std::cout<< feasible_range(p1,'<') <<std::endl;
        std::cout<< feasible_range(p1,'>') <<std::endl;
    }
    {
        std::map<clpoly::variable,clpoly::interval> tmp;
        tmp[y]=I1;tmp[z]=I1;
        clpoly::interval_poly<> p1=clpoly::assign(p,tmp);
        std::cout<<"p="<<p1<<std::endl;
        // {
        // clpoly::QQ l1=p1.back().second.get_l();
        // clpoly::QQ l2=p1.front().second.get_r();
        // std::cout<<"l1/l2="<<l1<<" "<<l2<<" "<<l1/l2<<std::endl;
        // } 
        std::cout<< feasible_range(p1,'=') <<std::endl;
        std::cout<< feasible_range(p1,'<') <<std::endl;
        std::cout<< feasible_range(p1,'>') <<std::endl;
    }
    
    return 0;
}
