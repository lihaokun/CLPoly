#include <iostream>
// #include <ginac/ginac.h>
// #include <time.h>
// using namespace std;
// using namespace GiNaC;
#include <boost/math/special_functions/prime.hpp>

int main()
{
    // symbol x("x");
    // symbol y("y");
    // symbol z("z");
    // symbol d("d");

    // // ex f=1+x+y+z+d;
    // // auto t=clock();
    // // ex f2=expand((f1+1)*(f1+2));
    // // std::cout<<"t:"<<double(clock()-t)/CLOCKS_PER_SEC << std::endl;
    // ex f = 4*x*y + x*z + 20*pow(y, 2) + 21*y*z + 4*pow(z, 2);
    // ex g = x*y + 3*x*z + 5*pow(y, 2) + 19*y*z + 12*pow(z, 2);
    // std::cout<<f<<std::endl;
    // std::cout<<g<<std::endl;    
    // std::cout<< resultant(f,g,x)<<std::endl;
    // std::cout<< expand((1+x+y+x*y)/(1+x))<<std::endl;
    std::cout<< boost::math::prime(0)<<std::endl;
    
    return 0;
}
