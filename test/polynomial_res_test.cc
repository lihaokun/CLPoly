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
    for (int i=0;i<100;++i)
    {
        //std::cout<<"test "<<i<<":\n";
        f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.02,10,-10);
        g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.02,10,-10);
        std::cout<<"f["<<i<<"]="<<f<<";"<<std::endl;
        std::cout<<"g["<<i<<"]="<<g<<";"<<std::endl; 
        t=clock();
        f1=resultant(f,g,x);
        s1+=double(clock()-t)/CLOCKS_PER_SEC;
        std::cout<<"o["<<i<<"]="<<f1<<";"<<std::endl; 
        std::cout<<"t["<<i<<"]="<<double(clock()-t)/CLOCKS_PER_SEC <<";"<<std::endl; 
        //std::cout<<"resultant    time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
        // t=clock();
        // f2=resultant_v1(f,g,x);
        // s2+=double(clock()-t)/CLOCKS_PER_SEC;
        // std::cout<<"resultant_v1 time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
        // if (f1==f2)
        //     std::cout<<"一致\n";
        // else
        // {
        //     std::cout<<"不一致\n";
        //     std::cout<<"f1="<<f1<<std::endl;
        //     std::cout<<"f2="<<f2<<std::endl;
             
             
        // }   
        // assert(f1==f2);
        
    }
    std::cout<<"(*resultant    total time:"<<s1<<"*)\n";
    // std::cout<<"resultant_v1 total time:"<<s2<<std::endl;
    return 0;
}
