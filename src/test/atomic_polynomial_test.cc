#include <iostream>
#include <vector>
#include "atomic_polynomial.hh"
#include <functional>
// template<typename Tm,typename Tc>
// std::ostream& operator<<(std::ostream& s, const clpoly::atomic_polynomial<Tm,Tc>& p) 
// {
//     s.put('[');
//     char comma[3] = {'\0', ' ', '\0'};
//     for (const auto& e : p) {
//         s << comma << e.first<<":"<<e.second;
//         comma[0] = ',';
//     }
//     return s << ']';
// }
int main(){
    clpoly::atomic_polynomial<int,int> p={{1,2},{3,4}};
    clpoly::atomic_polynomial<int,int> p2=p;
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p2:"<<p2<<std::endl;
    
    p.push_back({1,2});
    p.push_back({0,2});
    p2.push_back({2,5});
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p2:"<<p2<<std::endl;
    
    auto ptr=p.begin();
    std::cout<<"p.is_normal:"<<p.is_normal()<<std::endl;
    p.normalization();
    p2.normalization();
    std::cout<<"p2:"<<p2<<std::endl;
    std::cout<<"p:"<<p<<std::endl;
    clpoly::atomic_polynomial<int,int> p3(std::move(p));
    clpoly::atomic_polynomial<int,int> p1=std::move(p3);
    
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p1:"<<p1<<std::endl;
    std::cout<<"p3:"<<p3<<std::endl;
    
    std::cout<<"p1+p2:"<<(p1+p2)<<std::endl;
    std::cout<<"p2+p1:"<<(p2+p1)<<std::endl;
    std::cout<<"p1-p2:"<<(p1-p2)<<std::endl;
    std::cout<<"p2-p1:"<<(p2-p1)<<std::endl;
    p2.comp()=std::less<int>();
    std::cout<<"p2.is_normal:"<<p2.is_normal()<<std::endl;
    p2.normalization();
    std::cout<<"p2:"<<p2<<std::endl;
    p1.comp()=std::less<int>();
    p1.normalization();
    std::cout<<"p1+p2:"<<(p1+p2)<<"p1+p2.is_normal:"<<(p1+p2).is_normal()<<std::endl;

}