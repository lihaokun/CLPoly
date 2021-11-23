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
    clpoly::variable a("a");
    // for (size_t i=0;i<5;++i)
    // {
    //     auto f=clpoly::random_polynomial<clpoly::ZZ>({a,x,y,z},10,0.002,10,-10);
    //     std::cout<<f<<","<<std::endl;
    // }

    std::vector<clpoly::variable> l={x,y,z};
    std::vector<double> pl={0.5,0.1,0.02,0.01};
    auto p=pl.begin();
    for (auto i=l.begin();i!=l.end();++i)
    {
        auto j=i;
        ++j;
        assert(p!=pl.end());
        std::vector<clpoly::variable> l1(l.begin(),j);
        // std::cout<<l1<<","<<*p<<std::endl;
        auto f1=clpoly::random_polynomial<clpoly::ZZ>(l1,10,*p,10,-10);
        auto f=clpoly::random_polynomial<clpoly::ZZ>(l,10,0.005,10,-10);
        std::cout<<f*f1<<","<<std::endl;
        ++p;
    }
    return 0;
}