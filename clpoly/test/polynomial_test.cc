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
    clpoly::polynomial_ZZ p=2*y*z*x*x*z+d*d*d*y*z*z+1;
    std::cout<< p<<std::endl;
    auto l=p.variables();
    for (auto &i:l)
        std::cout<<i.first<<":"<<i.second<<" ";
    std::cout<<std::endl;
    clpoly::univariate_priority_order comp_z(z);
    clpoly::polynomial_<clpoly::ZZ,clpoly::univariate_priority_order> p2(&comp_z);
    clpoly::poly_convert(std::move(p),p2);
    std::cout<<p2<<std::endl;
    l=p2.variables();
    for (auto &i:l)
        std::cout<<i.first<<":"<<i.second<<" ";
    std::cout<<std::endl;
    //clpoly::__pair_vec_multiplies_compression_b=true;
    // std::cout<<pow(p2,2)<<std::endl;
    // std::cout<<"4*x^4*y^2*z^4+4*x^2*y^2*d^3*z^4+y^2*d^6*z^4+4*x^2*y*z^2+2*y*d^3*z^2+1\n";
    // std::cout<< p<<std::endl;
    // clpoly::polynomial_ZZ f=x*pow(y,2)+1;
    // clpoly::polynomial_ZZ g=2*pow(y,3)-pow(y,2)+pow(x,2)*y;

    clpoly::polynomial_ZZ f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},10,0.2,10,-10);
    clpoly::polynomial_ZZ g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},10,0.2,10,-10);
    // std::cout<<"g:="<<g<<":"<<std::endl;
    // std::cout<<"f:="<<f<<":"<<std::endl;
    auto t=clock();
    // std::cout<<"o:="<<clpoly::prem(g,f,y)<<":"<<std::endl;
    // std::cout<<"t:="<<(double(clock()-t)/CLOCKS_PER_SEC)<<";"<<std::endl;
    // std::cout<<"st := time():o1:= expand(prem(g, f, y)):time() - st;o1-o;"<<std::endl;
    
    // clpoly::variable r("r");
    // f=-pow(x,2)*pow(z,3) - pow(x,4) - pow(z,4) + pow(x,2) + 2*pow(z,2) - 1;
    // g=-pow(r,2)*pow(x,2) + pow(x,4) + pow(x,2)*pow(z,2) + pow(z,4) - 2*pow(z,2) + 1;
    g=-4*pow(x,10)+8*pow(x,9)*z-7*pow(x,8)*y*z-4*pow(x,7)*pow(z,3)+10*pow(x,4)*pow(y,2)*pow(z,4)-pow(x,4)*pow(z,6)-7*pow(x,3)*pow(y,6)*z+2*pow(x,3)*y*pow(z,6)-2*x*pow(y,4)*pow(z,5)-4*pow(y,8)*pow(z,2)-6*pow(y,7)*pow(z,3)+2*y*pow(z,9)-6*pow(x,7)*y*z+7*pow(x,5)*pow(y,3)*z-2*pow(x,3)*pow(y,2)*pow(z,4)-7*pow(x,2)*pow(y,5)*pow(z,2)-pow(x,2)*pow(y,2)*pow(z,5)-8*x*pow(y,4)*pow(z,4)-8*x*y*pow(z,7)-10*pow(y,8)*z+8*pow(x,8)+5*pow(x,7)*y+6*pow(x,7)*z-5*pow(x,5)*pow(y,3)-5*pow(x,4)*pow(y,4)+9*pow(x,3)*pow(y,3)*pow(z,2)-9*x*pow(y,6)*z+x*pow(y,4)*pow(z,3)-7*y*pow(z,7)-10*pow(x,5)*pow(z,2)+9*pow(x,3)*pow(y,4)+3*pow(x,3)*pow(y,2)*pow(z,2)+6*pow(x,2)*pow(y,3)*pow(z,2)-8*x*pow(y,4)*pow(z,2)+7*pow(y,2)*pow(z,5)-4*y*pow(z,6)+8*pow(x,4)*pow(y,2)-pow(x,4)*y*z-pow(x,2)*y*pow(z,3)+8*pow(x,2)*pow(z,4)-2*x*pow(y,5)-8*pow(y,6)+2*y*pow(z,5)+pow(x,2)*pow(y,2)*z+pow(x,2)*y*pow(z,2)+4*x*y*pow(z,3)+3*pow(y,4)*z-9*pow(x,4)+7*pow(x,3)*z+3*x*pow(z,3)+4*pow(y,2)*pow(z,2)-4*pow(z,4)-6*pow(x,2)*y+8*x*pow(y,2)+7*x*y*z-8*pow(y,3)-2*x*y-3*x*z-8*pow(z,2)+1;
    f=-8*pow(x,9)*z-pow(x,8)*pow(y,2)-2*pow(x,6)*pow(y,4)-6*pow(x,5)*pow(z,5)+2*pow(x,3)*pow(y,6)*z-3*pow(x,3)*pow(y,3)*pow(z,4)+2*pow(x,2)*pow(y,4)*pow(z,4)+9*pow(x,2)*pow(z,8)+3*x*pow(y,9)+3*x*pow(y,7)*pow(z,2)+5*pow(y,10)-5*pow(x,9)+5*pow(x,7)*pow(z,2)-6*pow(x,5)*pow(z,4)-3*pow(x,3)*pow(y,6)+pow(x,3)*pow(y,4)*pow(z,2)-2*pow(x,2)*pow(y,3)*pow(z,4)-4*pow(x,2)*pow(y,2)*pow(z,5)-8*pow(x,2)*y*pow(z,6)-3*x*pow(y,6)*pow(z,2)-9*x*pow(y,3)*pow(z,5)+2*pow(x,3)*pow(y,5)-pow(x,2)*pow(y,6)+4*x*pow(y,5)*pow(z,2)+2*pow(y,4)*pow(z,4)-3*pow(x,2)*pow(y,5)+pow(x,2)*pow(y,4)*z-8*pow(x,2)*pow(y,2)*pow(z,3)+8*pow(x,2)*pow(z,5)+5*x*pow(z,6)+7*pow(x,4)*pow(z,2)+5*pow(x,3)*y*pow(z,2)-7*pow(x,2)*pow(y,2)*pow(z,2)-5*x*pow(y,3)*pow(z,2)-4*x*y*pow(z,4)+3*pow(x,4)*y-3*pow(x,4)*z+pow(x,3)*y*z+6*pow(x,3)*pow(z,2)+9*pow(x,2)*pow(y,2)*z-2*x*pow(y,4)-9*x*pow(z,4)-6*pow(x,2)*y*z+5*pow(x,2)*pow(z,2)+9*x*y*pow(z,2)-8*pow(y,3)*z-pow(y,2)*pow(z,2)-9*x*pow(z,2)+y*z+z;

    // for (int i=0;i<10;++i)
    // {
    //     f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.2,10,-10);
    //     g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.2,10,-10);
        std::cout<<"g="<<g<<";"<<std::endl;
        std::cout<<"f="<<f<<";"<<std::endl;
        // f=-6*pow(x,4)-7*pow(x,2)*y-x*pow(y,3)*z-3*x*y-9;
        // g=10*pow(x,4)*y+2*pow(x,4)-7*pow(x,3)*pow(y,2)+9*pow(x,2)+9*pow(z,4)-9*y+2*z;
        // std::cout<<prem(f,g,x)<<std::endl;
        clpoly::__pair_vec_multiplies_compression_b=true;
        t=clock();
        std::cout<<"o="<<resultant(f,g,x)<<";"<<std::endl;
        std::cout<<"time:"<<(double(clock()-t)/CLOCKS_PER_SEC)<<std::endl;

        clpoly::__pair_vec_multiplies_compression_b=false;
        t=clock();
        std::cout<<"o="<<resultant(f,g,x)<<";"<<std::endl;
        std::cout<<"time:"<<(double(clock()-t)/CLOCKS_PER_SEC)<<std::endl;
    // }

    // clpoly::polynomial_ZZ p=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},5,0.05,10,-10);
    // std::cout<<"p="<<p<<std::endl;
    // auto l=p.variables();
    // for (auto &i:l)
    //     std::cout<<i.first<<":"<<i.second<<" ";
    // std::cout<<std::endl<<p.degree()<<std::endl;
    // t=clock();
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