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
    
    for (size_t i=0;i<5;++i)
    {
        auto f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.05,10,-10);
        std::cout<<f<<","<<std::endl;
    }
    return 0;
}