#include <iostream>
#include <vector>
#include <clpoly/clpoly.hh>
#include <functional>


int main(){
    clpoly::variable x("x");
    clpoly::variable x0("x0");
    clpoly::variable x1("x1");
    clpoly::variable x2("x2");
    clpoly::variable x3("x3");
    //clpoly::monomial px0=x0;
    //std::cout<<"x0:"<<px0<<std::endl;
    
    clpoly::monomial p={{x1,2},{x3,4}};
    std::vector<std::pair<clpoly::variable,int64_t>> v={{x0,1},{x1,2}};
    clpoly::monomial px1(std::move(v));
    std::cout<<"x0:"<<px1<<std::endl;
    
    std::cout<<"x0:"<<v.size()<<std::endl;
    
    clpoly::monomial p2=p;
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p2:"<<p2.str()<<std::endl;
    
    p.push_back({x1,2});
    p.push_back({x0,2});
    p2.push_back({x2,5});

    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p.deg:"<<p.deg()<<std::endl;
    
    std::cout<<"p2:"<<p2<<std::endl;
    
    auto ptr=p.begin();
    std::cout<<"p.is_normal:"<<p.is_normal()<<std::endl;
    p.normalization();
    p2.normalization();
    std::cout<<"p2:"<<p2<<std::endl;
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p.deg:"<<p.deg()<<std::endl;
    std::cout<<"p.deg(x1):"<<p.deg(x1)<<std::endl;
    clpoly::monomial p3(std::move(p));
    clpoly::monomial p1=std::move(p3);
    
    std::cout<<"p:"<<p<<" zore_check:"<<clpoly::zore_check<clpoly::monomial>()(p)<<std::endl;
    std::cout<<"p1:"<<p1<<std::endl;
    std::cout<<"p3:"<<p3<<std::endl;
    
    std::cout<<"p1*p2:"<<(p1*p2)<<std::endl;
    std::cout<<"p2*p1:"<<(p2*p1)<<std::endl;
    std::cout<<"p1/p2:"<<(p1/p2)<<std::endl;
    std::cout<<"p2/p1:"<<(p2/p1)<<std::endl;
    // clpoly::basic_monomial<clpoly::grevlex> m1=p1.data();
    // std::cout<<"m1:"<<m1<<std::endl;
    // clpoly::univariate_first_order x2f(x2);
    // std::cout<<"&x2f"<<&x2f<<std::endl;
    // clpoly::basic_monomial<clpoly::univariate_first_order> m1(&x2f);
    // m1=p1.data();
    // std::cout<<"m1:"<<m1<<std::endl;
    // clpoly::basic_monomial<clpoly::univariate_first_order> m2(&x2f);
    // m2=x2;
    // std::cout<<"m2:"<<m2<<std::endl;
    // clpoly::basic_monomial<clpoly::univariate_first_order> m3();
    // std::cout<<"m1*m2:"<<m1*m2<<std::endl;
    // std::cout<<"m1*m2?:"<<&((m1*m2).comp())<<std::endl;
    // std::cout<<"m1*m2?:"<<(m1*m2).is_normal()<<std::endl;
    
    //clpoly::basic_monomial<>
    // p2.comp(std::less<clpoly::variable>());
    // std::cout<<"p2.is_normal:"<<p2.is_normal()<<std::endl;
    // p2.normalization();
    // std::cout<<"p2:"<<p2<<std::endl;
    // p1.comp(std::greater<clpoly::variable>());
    // p1.normalization();
    // std::cout<<"p1*p2:"<<(p1*p2)<<" p1*p2.is_normal:"<<(p1*p2).is_normal()<<std::endl;
    // std::cout<<"p1.deg(x1)="<<p1.deg(x1)<<std::endl;
    
}