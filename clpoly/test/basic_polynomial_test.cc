#include <iostream>
#include <vector>
#include "basic_polynomial.hh"
#include "monomial.hh"
#include "variable.hh"
#include "polynomial.hh"
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
    clpoly::monomial m1={{x,1}};
    std::cout<<clpoly::_monomial_compression(m1,{x,y,z,g})<<std::endl;
    clpoly::_monomial_decompression(clpoly::_monomial_compression(m1,{x,y,z,g}),m1,{x,y,z,g},m1.comp_ptr());
    std::cout<<m1.deg()<<std::endl;
    std::cout<<"p.power(2):"<<pow(p,2)<<std::endl;
    std::cout<<"p.power(2):x^2+2*x*y+2*x*z+2*x*g+y^2+2*y*z+2*y*g+z^2+2*z*g+g^2+2*x+2*y+2*z+2*g+1\n";
    std::cout<<"(p+1)*(p+2):"<<((p+n1)*(p+n2))<<std::endl;
    std::cout<<"(p+1)*(p+2):x^2+2*x*y+2*x*z+2*x*g+y^2+2*y*z+2*y*g+z^2+2*z*g+g^2+5*x+5*y+5*z+5*g+6\n";
    //std::cout<<"p.power(2).coef(x*y):"<<pow(p,2).coef({{x,1},{y,1}}))<<std::endl;
    std::cout<<"(p+1)*(p+2)/(p+1):"<<((p+n1)*(p+n2))/(p+n1)<<std::endl;
    // clpoly::basic_polynomial<clpoly::monomial,clpoly::QQ> p4,p5;
    // //clpoly::pair_vec_div(p4.data(),((p+n1)*(p+n2)).data(),(p*2+n1).data(),(p+n1).comp());
    // clpoly::polynomial_div(p4,(p+1)*(p+2),(p*2+1));
    // std::cout<<"(p+1)*(p+2)/(2*p+1):"<<p4<<std::endl;
    // clpoly::poly_convert((p*2+1),p5);
    // std::cout<<p5<<std::endl;
    // std::cout<<p4*p5<<std::endl;
    
    // clpoly::variable x1("x1");
    // clpoly::variable x2("x2");
    // clpoly::variable x3("x3");
    // clpoly::variable x4("x4");
    // clpoly::variable x5("x5");
    // clpoly::variable x6("x6");
    // clpoly::polynomial_ZZ f1;
    // clpoly::polynomial_ZZ f2;
    
    // double t1=0,t2=0;
    // auto t=clock();
    // for(int i=0;i<1;++i)
    // {
    //     f1=clpoly::random_polynomial<clpoly::ZZ>({x1},500,0.01,10,-10);
    //     f2=clpoly::random_polynomial<clpoly::ZZ>({x1},500,0.01,10,-10);
    //     std::cout<<"f1.size:"<<f1.size()<<std::endl;
    //     std::cout<<"f2.size:"<<f2.size()<<std::endl;
    //     clpoly::__pair_vec_multiplies_compression_b=true;
    //     t=clock();
    //     std::cout<<(f1*f2)<<std::endl;
    //     t1+=((double)clock()-t)/CLOCKS_PER_SEC;
    //     printf ("(%f seconds).\n",((double)clock()-t)/CLOCKS_PER_SEC);
    //     clpoly::__pair_vec_multiplies_compression_b=false;
    //     t=clock();
    //     std::cout<<(f1*f2)<<std::endl;
    //     t2+=((double)clock()-t)/CLOCKS_PER_SEC;
    //     printf ("(%f seconds).\n",((double)clock()-t)/CLOCKS_PER_SEC);
    // }
    // //std::cout<<t1<<"\n"<<t2<<"\n"<<t2-t1<<std::endl;
    // clpoly::__pair_vec_multiplies_compression_b=true;

    clpoly::basic_polynomial<clpoly::monomial,clpoly::ZZ> p2=pow(p,20);
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
    t=clock();
    std::cout<<((p2+n1)==(p3/(p2+n2)))<<std::endl;
    printf ("(%f seconds).\n",((float)clock()-t)/CLOCKS_PER_SEC);
    

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