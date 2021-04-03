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
    clpoly::polynomial_ZZ f1,f2,f3,f4,f,g,G,r;
    // f=4*pow(z,10)+6*pow(x,9)+6*pow(y,2)*pow(z,5)+2*pow(y,4)-5*x-6;
    // g=5*x*pow(y,5)*pow(z,4)+9*pow(y,2)*pow(z,8)-5*pow(x,2)*pow(y,5)*pow(z,2)+6*x*y*pow(z,4)-5;
    // std::cout<<"f="<<f<<std::endl;
    // std::cout<<"g="<<g<<std::endl;
    // // // // std::cout<<"prem(f,g,x)="<<prem(f,g,x)<<std::endl;
    // std::cout<<"prem_v2(f,g,x)="<<prem_v2(f,g,x)<<std::endl;
    
    
    time_t t;
    double s1=0,s2=0;
    int n;
    std::string s;
    std::cout<<"Test round:";
    std::cin>>s;
    n=std::stoi(s); 
    for (int i=0;i<n;++i)
    {
        std::cout<<"test "<<i<<":\n";
        f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.01,10,-10);
        g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.01,10,-10);
        std::cout<<"f="<<f<<";"<<std::endl;
        std::cout<<"g="<<g<<";"<<std::endl; 
        if (g.empty())
            continue;
        t=clock();
        // f1=clpoly::pquo(r,f,g,x);
        f1=clpoly::pquo(f,g,x);
        // std::cout<<"q="<<f1<<";"<<std::endl; 
        // std::cout<<"r="<<r<<";"<<std::endl; 
        s1+=double(clock()-t)/CLOCKS_PER_SEC;
        // std::cout<<"pquo    time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
        t=clock();
        
        f2=clpoly::prem(f,g,x);
        // std::cout<<"pquo    time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
        // std::cout<<(f2==r)<<std::endl;
        auto d=clpoly::degree(f,x)-clpoly::degree(g,x)+1;
        if (d<=0)
        {
            d=0;
        }
        std::cout<<(f*pow(clpoly::leadcoeff(g,x),d)==(f2+f1*g))<<std::endl;
        
        
        
        // t=clock();
        // f2=prem_v1(f,g,x);
        // s2+=double(clock()-t)/CLOCKS_PER_SEC;
        // std::cout<<"prem_v1 time:"<<double(clock()-t)/CLOCKS_PER_SEC<<std::endl;
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
    std::cout<<"pquo    total time:"<<s1<<std::endl;
    // std::cout<<"prem_v1 total time:"<<s2<<std::endl;
    
    return 0;
}
