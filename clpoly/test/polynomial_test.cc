#include "polynomial.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
clpoly::polynomial_ZZ read_file(std::string s)
{
    long numn,numv;
    std::ifstream fin(s);
    fin>>numn>>numv;
    std::vector<clpoly::variable> var;
    var.reserve(numv);
    for (char i='a';i<'a'+numv;var.push_back(std::string(1,i)),++i);
    std::vector<std::pair<clpoly::variable,int64_t>> m;
    std::vector<std::pair<clpoly::monomial,clpoly::ZZ>> p;
    p.reserve(numn);
    long tmp;
    for (long i=0;i<numn;++i)
    {
        m.reserve(numv);
        for (auto &j:var)
        {
            fin>>tmp;
            if (tmp!=0)
                m.push_back(std::pair<clpoly::variable,int64_t>(j,tmp));
        }
        fin>>tmp;
        if (tmp!=0)
            p.push_back(std::pair<clpoly::monomial,clpoly::ZZ>(std::move(m),clpoly::ZZ(tmp)));
    }
    fin.close();
    return clpoly::polynomial_ZZ(p);
}
int main(){
    clpoly::variable x("x");
    clpoly::variable y("y");
    clpoly::variable z("z");
    clpoly::variable d("d");
    clpoly::polynomial_ZZ p=2*x*z*x*x+d*y*z*z+1;
    std::cout<< p<<std::endl;
    clpoly::univariate_first_order comp_z(z);
    clpoly::polynomial_<clpoly::ZZ,clpoly::univariate_first_order> p2(&comp_z);
    clpoly::poly_convert(std::move(p),p2);
    std::cout<<p2<<std::endl;
    std::cout<< p<<std::endl;
    clpoly::polynomial_ZZ f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},40,0.2,10,-10);
    clpoly::polynomial_ZZ g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},40,0.2,10,-10);
    std::cout<<"f:="<<f<<";"<<std::endl;
    std::cout<<"g:="<<g<<";"<<std::endl;
    auto t=clock();
    std::cout<<"o1:="<<clpoly::prem(g,f,z)<<";"<<std::endl;
    std::cout<<"time:"<<(double(clock()-t)/CLOCKS_PER_SEC)<<std::endl;
    //std::cout<<clpoly::random_polynomial<clpoly::ZZ>({x,y,z},5,0.2,10,-10)<<std::endl;
    // auto l=p.variables();
    // for (auto &i:l)
    //     std::cout<<i.first<<":"<<i.second<<" ";
    // std::cout<<std::endl<<p.degree()<<std::endl;
    // auto t=clock();
    // auto PP=read_file("/home/ker/Documents/Bigpoly/j621_data.txt");
    // std::cout<<PP.size()<<std::endl;
    // std::cout<<"( "<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";
    // t=clock();
    // auto ll=PP.variables();
    // for (auto &i:ll)
    //     std::cout<<i.first<<":"<<i.second<<" ";
    // std::cout<<std::endl;
    // std::cout<<"("<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";
    // t=clock();
    // std::cout<<PP.degree()<<std::endl;
    // std::cout<<"("<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";


    return 0;
}