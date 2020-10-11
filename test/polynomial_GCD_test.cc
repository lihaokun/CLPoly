#include <clpoly.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
int main(int argc, char const *argv[])
{
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::polynomial_ZZ f1,f2,f3,f4,f,g,G;

    
    time_t t;
    double s1=0,s2=0;
    std::string s;
    std::cout<<"Test round:";
    std::cin>>s;
    int n=std::stoi(s);
    for (int i=0;i<n;++i)
    {
        //std::cout<<"test "<<i<<":\n";
        f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.5,10,-10);
        g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.5,10,-10);
        std::cout<<"f["<<i+1<<"]="<<f<<";"<<std::endl;
        std::cout<<"g["<<i+1<<"]="<<g<<";"<<std::endl; 
        t=clock();
        f1=clpoly::polynomial_GCD(f*f,g*f);
        s1+=double(clock()-t)/CLOCKS_PER_SEC;
        std::cout<<"o["<<i+1<<"]="<<f1<<";"<<std::endl; 
        std::cout<<"t["<<i+1<<"]="<<double(clock()-t)/CLOCKS_PER_SEC <<";"<<std::endl; 
        
    }
    std::cout<<"(*GCD    total time:"<<s1<<"*)\n";
    // std::cout<<"resultant_v1 total time:"<<s2<<std::endl;
    return 0;
}