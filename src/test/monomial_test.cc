#include <iostream>
#include <vector>
#include "monomial.hh"
#include "variable.hh"
#include <functional>

std::ostream& operator<<(std::ostream& s, const clpoly::monomial& p) 
{
    s.put('[');
    char comma[3] = {'\0', ' ', '\0'};
    for (const auto& e : p) {
        s << comma << e.first<<":"<<e.second;
        comma[0] = ',';
    }
    return s << ']';
}
int main(){
    clpoly::variable x0("x0");
    clpoly::variable x1("x1");
    clpoly::variable x2("x2");
    clpoly::variable x3("x3");
    
    clpoly::monomial p={{x1,2},{x3,4}};
    clpoly::monomial p2=p;
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p2:"<<p2<<std::endl;
    
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
    
    clpoly::monomial p3(std::move(p));
    clpoly::monomial p1=std::move(p3);
    
    std::cout<<"p:"<<p<<std::endl;
    std::cout<<"p1:"<<p1<<std::endl;
    std::cout<<"p3:"<<p3<<std::endl;
    
    std::cout<<"p1+p2:"<<(p1+p2)<<std::endl;
    std::cout<<"p2+p1:"<<(p2+p1)<<std::endl;
    std::cout<<"p1-p2:"<<(p1-p2)<<std::endl;
    std::cout<<"p2-p1:"<<(p2-p1)<<std::endl;
    p2.comp()=std::greater<clpoly::variable>();
    std::cout<<"p2.is_normal:"<<p2.is_normal()<<std::endl;
    p2.normalization();
    std::cout<<"p2:"<<p2<<std::endl;
    p1.comp()=std::greater<clpoly::variable>();
    p1.normalization();
    std::cout<<"p1+p2:"<<(p1+p2)<<"p1+p2.is_normal:"<<(p1+p2).is_normal()<<std::endl;

}