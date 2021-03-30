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
    clpoly::variable d("d");
        clpoly::variable r("r");
    clpoly::polynomial_ZZ f1,f2,f3,f4,f,g,G;

    // f=-3*pow(x,4)*pow(y,6)-10*pow(x,9)-2*pow(x,4)*pow(y,2)+9;
    // g=5*pow(x,9)+pow(x,2)*pow(y,6)-10;
        
    time_t t;
    double s1=0,s2=0;
    for (int i=0;i<100;++i)
    {
        f=clpoly::random_polynomial<clpoly::ZZ>({x,y},10,0.1,10,-10);
        g=clpoly::random_polynomial<clpoly::ZZ>({x,y},10,0.1,10,-10);
        if (degree(f,x)<degree(g,x))
            std::swap(f,g);
        if (degree(g,x)==0)
            continue;
        std::cout<<"f["<<i+1<<"]="<<f<<";"<<std::endl;
        std::cout<<"g["<<i+1<<"]="<<g<<";"<<std::endl; 
        t=clock();
        auto o=subresultant(f,g,x);
        s1+=double(clock()-t)/CLOCKS_PER_SEC;
        std::cout<<"o["<<i+1<<"]=List"<<o<<";"<<std::endl; 
        std::cout<<"t["<<i+1<<"]="<<double(clock()-t)/CLOCKS_PER_SEC <<";"<<std::endl; 
        
    }
    std::cout<<"(*resultant    total time:"<<s1<<"*)\n";
    return 0;


}
