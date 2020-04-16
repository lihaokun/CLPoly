#include <iostream>
#include <vector>
#include "basic_polynomial.hh"
#include "monomial.hh"
#include "variable.hh"
#include "number.hh"
#include <functional>
#include <time.h>

int main(){
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::variable g("g");
    //clpoly::monomial m1={{x,1}};
    clpoly::basic_polynomial<clpoly::monomial,clpoly::ZZ> p={{{{x,1}},1},{{{y,1}},1},{{{z,1}},1},{{{g,1}},1},{{},1}};
    std::cout<<p<<" normal?"<<p.is_normal()<<std::endl;
    
    p.normalization();
    std::cout<<p<<std::endl;
    clpoly::basic_polynomial<clpoly::monomial,clpoly::ZZ> n1={{{},1}};
    clpoly::basic_polynomial<clpoly::monomial,clpoly::ZZ> n2={{{},2}};
    //std::cout<<"n1<n2?"<<(n1.begin()->first<n2.begin()->first)<<std::endl;
    // // clpoly::basic_polynomial<int,int> p={{1,2},{3,4}};
    // clpoly::basic_polynomial<int,int> p2=p;
    std::cout<<"p.power(2):"<<(p.power(2))<<std::endl;
    std::cout<<"(p+1)*(p+2):"<<((p+n1)*(p+n2))<<std::endl;
    std::cout<<"p.power(2).coef(x*y):"<<(p.power(2).coef({{x,1},{y,1}}))<<std::endl;
    clpoly::basic_polynomial<clpoly::monomial,clpoly::ZZ> p2=p.power(20);
    auto t=clock();
    auto p3=(p2+n1)*(p2+n2);
    printf ("(%f seconds).\n",((float)clock()-t)/CLOCKS_PER_SEC);
    std::cout<<112911876*1.0/p3.size()<<std::endl;//831.757
    auto k=p3.begin()->second;
    for(auto & i:p3)
        if (k<i.second)
            k=i.second;
    std::cout<<k<<std::endl;
    std::cout<<"7656714453153197981835000"<<std::endl;
    // p.push_back({1,2});
    // p.push_back({0,2});
    // p2.push_back({2,5});
    // std::cout<<"p:"<<p<<std::endl;
    // std::cout<<"p2:"<<p2<<std::endl;
    
    // auto ptr=p.begin();
    // std::cout<<"p.is_normal:"<<p.is_normal()<<std::endl;
    // p.normalization();
    // p2.normalization();
    // std::cout<<"p2:"<<p2<<std::endl;
    // std::cout<<"p:"<<p<<std::endl;
    // clpoly::basic_polynomial<int,int> p3(std::move(p));
    // clpoly::basic_polynomial<int,int> p1=std::move(p3);
    
    // std::cout<<"p:"<<p<<std::endl;
    // std::cout<<"p1:"<<p1<<std::endl;
    // std::cout<<"p3:"<<p3<<std::endl;
    
    // std::cout<<"p1+p2:"<<(p1+p2)<<std::endl;
    // std::cout<<"p2+p1:"<<(p2+p1)<<std::endl;
    // std::cout<<"p1-p2:"<<(p1-p2)<<std::endl;
    // std::cout<<"p2-p1:"<<(p2-p1)<<std::endl;
    // p2.comp()=std::less<int>();
    // std::cout<<"p2.is_normal:"<<p2.is_normal()<<std::endl;
    // p2.normalization();
    // std::cout<<"p2:"<<p2<<std::endl;
    // p1.comp()=std::less<int>();
    // p1.normalization();
    // std::cout<<"p1+p2:"<<(p1+p2)<<"p1+p2.is_normal:"<<(p1+p2).is_normal()<<std::endl;

}