#include <clpoly/variable.hh>
#include <iostream>
#include <map>
#include <unordered_map>

int main(){
    clpoly::variable v;
    std::cout<<"v "<<v.serial()<<" "<<v<<std::endl;
    v=clpoly::variable::new_variable("x");
    std::cout<<"x1 "<<clpoly::variable::get_variable("x1").serial()<<std::endl;
    clpoly::variable v1=v;
    std::cout<<"v1 "<<v1.serial()<<" "<<v1<<std::endl;
    clpoly::variable v2("A");
    clpoly::variable v3=clpoly::variable::get_variable("A");
    std::cout<<"v2 "<<v2.serial()<<" "<<v2<<std::endl;
    //std::cout<<"v2 "<<std::size_t(v2)<<" "<<v2<<std::endl;
    
    std::cout<<"v3 "<<v3.serial()<<" "<<v3<<std::endl;
    clpoly::variable::del_variable("x1");
    
    clpoly::variable v4=clpoly::variable::get_variable(3);
    std::cout<<"v4 "<<v4.serial()<<" "<<v4<<std::endl;
    std::cout<<"v!=v1 "<<(v!=v1)<<" v<v4 "<<(v<v4)<<std::endl;
    v=clpoly::variable::new_variable();
    std::cout<<"v "<<v.serial()<<" "<<v<<std::endl;
    std::map<clpoly::variable,int> dct={
    {v1,1},{v4,2}
    };
    for (auto & i:dct)
    {
        std::cout<<i.first<<":"<<i.second<<std::endl;
    }
    std::unordered_map<clpoly::variable,int> dct1={
    {v1,1},{v4,2}
    };
    for (auto & i:dct1)
    {
        std::cout<<i.first<<":"<<i.second<<std::endl;
    }

}