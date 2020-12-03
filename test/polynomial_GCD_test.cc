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
    clpoly::polynomial_ZZ f2,f3,f4,f,g,G;

    // f=pow(x,4)+25*pow(x,3)+145*pow(x,2)-171*x-360;
    // g=pow(x,5)+14*pow(x,4)+15*pow(x,3)-pow(x,2)-14*x-15;
    // std::cout<<clpoly::polynomial_GCD(f,g)<<std::endl;
    // f=2*pow(x,4)-7*pow(x,3)-4*pow(x,2)-4*x-15;
    // g=4*pow(x,5)+4*pow(x,3)-7*pow(x,2)-2*pow(x,4)+x-12;
    // std::cout<<clpoly::polynomial_GCD(f,g)<<std::endl;
    // f=8*pow(x,2)*y*pow(z,2)-8*pow(x,2)*pow(z,3)+3*x*pow(y,3)*z+4*x*y*pow(z,2)+5*pow(x,2)*z-9*x*pow(y,2)-6*pow(y,3);
    // g=2*pow(x,5)-6*x*pow(y,4)-2*pow(x,2)*y*z-4*pow(x,2)*z+3*pow(z,2)-1;
    // std::cout<<clpoly::polynomial_GCD(f*f,g*f)<<std::endl;
    // f=-pow(x,2)*pow(y,2)*z+7*x*y*pow(z,2)-9*x*pow(z,2)-3*pow(y,3)+4*pow(y,2)+5;
    // g=-6*pow(x,3)*y*z-8*pow(z,3)+5*z-10;
    // std::cout<<clpoly::polynomial_GCD(f*f,g*f)<<std::endl;
    // f=4*pow(x,4)*y-3*pow(x,3)*y*z-3*pow(x,2)*pow(y,3)+x*pow(y,4)+5*y*pow(z,4)-10*pow(z,5)-7*x*pow(y,2)*z-x*y*pow(z,2)-pow(x,2)*z-10*x*y*z+3*y*pow(z,2)-9*pow(x,2)-3;
    // g=-6*pow(x,2)*pow(y,2)*z+2*x*pow(y,3)*z-3*x*y*pow(z,3)-3*x*y*pow(z,2)+8*y*pow(z,3)+10*x*pow(y,2)+x*pow(z,2)+10*y*pow(z,2)-7;
    // std::cout<<clpoly::polynomial_GCD(f*f,g*f)<<std::endl;
    // f=2*y;
    // g=10*pow(y,2)*pow(z,4)+6*pow(y,2)*pow(z,2)-20*y*pow(z,5)-6*y;
    // std::cout<<clpoly::polynomial_GCD(f,g)<<std::endl;
    f=pow(x,5)-3*pow(x,4)+4*pow(x,3)-4*pow(x,2)+3*x-1;
    clpoly::lex_<clpoly::custom_var_order> mo(clpoly::custom_var_order({z,y,x}));
    clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> f_(&mo); 
    clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> g_(&mo);
    
    clpoly::poly_convert(f,f_);
    auto l=clpoly::squarefree(f);
    std::cout<<f<<":";
    for(auto &i:l)
    {
        std::cout<<"{"<<i.first<<","<<i.second<<"} ";
    }
    std::cout<<std::endl;

    auto l1=clpoly::squarefreebasis(std::vector<clpoly::polynomial_ZZ>({f,f}));
    std::cout<<f<<":";
    for(auto &i:l1)
    {
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
    
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
        clpoly::poly_convert(f,f_);
        clpoly::poly_convert(g,g_);
        std::cout<<"f["<<i+1<<"]="<<f_<<";"<<std::endl;
        std::cout<<"g["<<i+1<<"]="<<g_<<";"<<std::endl; 
        t=clock();
        auto f1=clpoly::polynomial_GCD(f_*f_,g_*f_);
        //  auto f1=clpoly::polynomial_GCD(f*f,g*f);
        s1+=double(clock()-t)/CLOCKS_PER_SEC;
        std::cout<<"o["<<i+1<<"]="<<f1<<";"<<std::endl; 
        std::cout<<"t["<<i+1<<"]="<<double(clock()-t)/CLOCKS_PER_SEC <<";"<<std::endl; 
        
    }
    std::cout<<"(*GCD    total time:"<<s1<<"*)\n";
    // std::cout<<"resultant_v1 total time:"<<s2<<std::endl;
    return 0;
}