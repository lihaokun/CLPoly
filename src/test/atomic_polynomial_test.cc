#include <iostream>
#include <vector>
#include "atomic_polynomial.hh"
template<typename Tm,typename Tc>
std::ostream& operator<<(std::ostream& s, const clpoly::atomic_polynomial<Tm,Tc>& p) 
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
    clpoly::atomic_polynomial<int,int> p={{1,2},{3,4}};
    p.push_back({1,2});
    std::cout<<p;

}