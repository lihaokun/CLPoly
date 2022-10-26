#include <clpoly/clpoly.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
int main(int argc, char const *argv[])
{
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::variable a1("a1");
    clpoly::variable a2("a2");
    clpoly::variable a3("a3");
    clpoly::variable a4("a4");
    
    std::vector<clpoly::variable> l={x,y,z,a1};
    // std::vector<double> pl={0.5,0.1,0.02,0.01};
    // auto p=pl.begin();
    // for (auto i=0;i!=4;++i)
    // {
    //     // auto j=i;
    //     // ++j;
    //     // assert(p!=pl.end());
    //     // std::vector<clpoly::variable> l1(l.begin(),j);
    //     // std::cout<<l1<<","<<*p<<std::endl;
    //     // auto f1=clpoly::random_polynomial<clpoly::ZZ>(l1,10,*p,10,-10);
    //     auto f=clpoly::random_polynomial<clpoly::ZZ>(l,5,0.04,10,-10);
    //     std::cout<<f<<","<<std::endl;
    //     // ++p;
    // }
    auto f=clpoly::random_polynomial<clpoly::ZZ>({a1,x},10,0.1,10,-10);
    std::cout<<f<<","<<std::endl;
    f=clpoly::random_polynomial<clpoly::ZZ>({a1,x,y},10,0.01,10,-10);
    std::cout<<f<<","<<std::endl;
    // auto f=clpoly::random_polynomial<clpoly::ZZ>({},10,0.001,10,-10);
    // std::cout<<f<<","<<std::endl;
    
    return 0;
}