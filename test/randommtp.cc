#include <clpoly/clpoly.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
int main(int argc, char const *argv[])
{
    clpoly::variable x("x");
    clpoly::variable y("sin(x)");
    clpoly::variable z("cos(x)");
    // for (size_t i=0;i<5;++i)
    // {
    //     auto f=clpoly::random_polynomial<clpoly::ZZ>({a,x,y,z},10,0.002,10,-10);
    //     std::cout<<f<<","<<std::endl;
    // }

    std::vector<clpoly::variable> l={x,y,z};
    std::vector<std::pair<int,double>> pl={{5,0.1},{10,0.015},{10,0.05},{10,0.1},{15,0.01},{15,0.05},{15,0.1},
    {20,0.0025}};
    // auto p=pl.begin();
    int n=1000;
    for (auto p:pl)
    {
        std::cout<<"# "<<p.first<<" "<<p.second<<" "<<n<<std::endl;
        for (auto i=0;i!=10;++i)
        {
            // auto j=i;
            // ++j;
            // assert(p!=pl.end());
            // std::vector<clpoly::variable> l1(l.begin(),j);
            // std::cout<<l1<<","<<*p<<std::endl;
            // auto f1=clpoly::random_polynomial<clpoly::ZZ>(l1,10,*p,10,-10);
            auto f=clpoly::random_polynomial<clpoly::ZZ>(l,p.first,p.second,n,-n);
            std::cout<<"[ "<<f<<" , "<<x<<" ], "<<std::endl;
            // ++p;
        }
    }
    return 0;
}