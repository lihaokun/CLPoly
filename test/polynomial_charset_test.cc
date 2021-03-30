#include <clpoly/clpoly.hh>
#include <iostream>

int main(int argc, char const *argv[])
{
    
       clpoly::variable x1("x1"),x2("x2"),x3("x3"),x4("x4");
    clpoly::lex_<clpoly::custom_var_order> mo(clpoly::custom_var_order({x4,x3,x2,x1}));
    clpoly::polynomial_ZZ p1,p2,p3,p4;
    p1=x1*pow(x4,2)+pow(x4,2)-x1*x2*x4-x2*x4+x1*x2+3*x2;
    p2=x1*x4+x3-x1*x2;
    p3=x3*x4-2*pow(x2,2)-x1*x2-1;
    
    
    clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> p1_(&mo),p2_(&mo),p3_(&mo);
    std::vector<clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>>> P;
    clpoly::poly_convert(p1,p1_);
    clpoly::poly_convert(p2,p2_);
    clpoly::poly_convert(p3,p3_);
    auto p4_=prem(p1_,p2_,x4);
    auto p5_=prem(p3_,p2_,x4);
    //std::cout<<"p1_="<<p1_<<";"<<std::endl;
    // std::cout<<"p2_="<<p2_<<";"<<std::endl;
    // std::cout<<"p3_="<<p3_<<";"<<std::endl;
    // std::cout<<"p4_="<<p4_<<";"<<std::endl;
    //std::cout<<"p5_="<<p5_<<";"<<std::endl;
    //std::cout<<"fist_var="<<(clpoly::get_first_var(p1_))<<";"<<std::endl;
    // std::cout<<"rank(p2_)<rank(p1_)?"<<rankcompare(p2_,p1_)<<";"<<std::endl;
    // std::cout<<rankcompare(p3_,p2_)<<";"<<rankcompare(p2_,p3_)<<std::endl;
    //std::cout<<"p3_ reduce p1_?"<<is_reduced(p3_,p1_)<<std::endl;
    //std::cout<<"p1_ reduce p3_?"<<is_reduced(p1_,p3_)<<std::endl;
    P.push_back(p1_);
    P.push_back(p2_);
    P.push_back(p3_);
    // P.push_back(p4_);
    // P.push_back(p5_);
    auto C=clpoly::charset(P); 
    std::cout<<C<<std::endl;
    return 0;
    return 0;
}