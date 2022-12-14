#include <clpoly/clpoly.hh>
#include <boost/container_hash/hash.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
clpoly::polynomial_ZZ read_file(std::string s)
{
    long numn,numv;
    std::ifstream fin(s);
    if (!fin)
    {
        // std::cout<<"无法读取文件!\n";
        throw "无法读取文件!";
    }
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
    clpoly::variable r("r");
    clpoly::variable x1("x1"),x2("x2"),x3("x3");
    clpoly::variable x7("x7"),x8("x8");
    clpoly::polynomial_ZZ p={{{{x,1}},1}};
    clpoly::lex_<clpoly::custom_var_order> mo;
    mo=clpoly::lex_<clpoly::custom_var_order>(clpoly::custom_var_order({x8,x7,x1}));
    clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> p_1(&mo);
    std::cout<<"p="<< p<<std::endl;
    p=1;
    p=2*(2-x1)*x8+x7-2;
    std::cout<<"p="<< p<<std::endl;
    p_1=p;
    std::cout<<"p_1="<< p_1<<std::endl;
    
    clpoly::monomial m=pow(x,3);
    p=z*z;
    time_t t;
    std::cout<< p<<std::endl;
    
    std::cout<< 2*x*x*y+1<<std::endl;
    // auto l=p.variables();
    // for (auto &i:l)
    //     std::cout<<i.first<<":"<<i.second<<" ";
    // std::cout<<std::endl;
    // clpoly::univariate_priority_order comp_z(z);
    // clpoly::polynomial_<clpoly::ZZ,clpoly::univariate_priority_order> p2(&comp_z);
    // clpoly::poly_convert(std::move(p),p2);
    // std::cout<<p2<<std::endl;
    // l=p2.variables();
    // for (auto &i:l)
    //     std::cout<<i.first<<":"<<i.second<<" ";
    // std::cout<<std::endl;
    //clpoly::__is_monomial_compression=true;
    // std::cout<<pow(p2,2)<<std::endl;
    // std::cout<<"4*pow(x,4)*pow(y,2)*pow(z,4)+4*pow(x,2)*pow(y,2)*d^3*pow(z,4)+pow(y,2)*d^6*pow(z,4)+4*pow(x,2)*y*pow(z,2)+2*y*d^3*pow(z,2)+1\n";
    // std::cout<< p<<std::endl;
    clpoly::polynomial_ZZ f=x*pow(y,2);
    clpoly::polynomial_ZZ g=2*pow(y,3)-pow(y,2)+pow(x,2)*y;
    
    f=2*(2-7*x1+pow(x1,2)*x2)-(x3-x1);
    std::cout<<"f:="<<f<<":"<<std::endl;
    f=-pow(x,2)*pow(z,3) - pow(x,4) - pow(z,4) + pow(x,2) + 2*pow(z,2) - 1;
    g=-pow(r,2)*pow(x,2) + pow(x,4) + pow(x,2)*pow(z,2) + pow(z,4) +x- 2*pow(z,2) + 1;
    // clpoly::polynomial_ZZ f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},10,0.2,10,-10);
    // clpoly::polynomial_ZZ g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},10,0.2,10,-10);
    clpoly::polynomial_ZZ o;
     clpoly::polynomial_ZZ o1;
    std::cout<<"f:="<<f<<":"<<std::endl;
    std::cout<<"g:="<<g<<":"<<std::endl;
    std::cout<<clpoly::coeff(g,z)<<std::endl;
    // std::cout<<"o:="<<clpoly::prem(g,f,y)<<":"<<std::endl;
    // std::cout<<"t:="<<(double(clock()-t)/CLOCKS_PER_SEC)<<";"<<std::endl;
    // std::cout<<"st := time():o1:= expand(prem(g, f, y)):time() - st;o1-o;"<<std::endl;
    
 

    // std::cout<<"t=Association[];t2=Association[];\n";


    // clpoly::polynomial_ZZ p=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},5,0.05,10,-10);
    // std::cout<<"p="<<p<<std::endl;
    // auto l=p.variables();
    // for (auto &i:l)
    //     std::cout<<i.first<<":"<<i.second<<" ";
    // std::cout<<std::endl<<p.degree()<<std::endl;
    t=clock();
    clpoly::polynomial_ZZ PP;
    try
    {   
        PP=read_file("j621_data.txt");
    }
    catch(const char* msg)
    {
        std::cout<<msg<<std::endl;
        return 1;
    }
    std::cout<<PP.size()<<std::endl;
    std::cout<<"( "<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";
    t=clock();
    auto ll=PP.variables();
    for (auto &i:ll)
        std::cout<<i.first<<":"<<i.second<<" ";
    std::cout<<std::endl;
    std::cout<<"("<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";
     t=clock();
    std::cout<<std::hash<clpoly::polynomial_ZZ>()(PP)<<std::endl;
    std::cout<<"("<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";
    // t=clock();
    // std::cout<<PP.degree()<<std::endl;
    // std::cout<<"("<<double(clock()-t)/CLOCKS_PER_SEC<<"s)\n";

    // clpoly::Zp a(10,7);
    // clpoly::Zp b(-10,7);
    // a=-2;
    // std::cout<<a<<" "<<b<<" "<<a*b<<" "<<a/b<<" "<<-a<<std::endl;

    // g=-pow(r,2)*pow(x,2) + pow(x,4) + pow(x,2)*pow(z,2) + pow(z,4) - 2*pow(z,2) + 1;
    // std::cout<<g<<std::endl;
    // std::cout<<clpoly::polynomial_mod(std::move(g),7)<<std::endl;
    // a=std::move(b);
    // std::cout<<a<<b<<std::endl;
    // f=pow(x,4)+25*pow(x,3)+145*pow(x,2)-171*x-360;
    // g=pow(x,5)+14*pow(x,4)+15*pow(x,3)-pow(x,2)-14*x-15;
    // f=2*pow(x,4)-7*pow(x,3)-4*pow(x,2)-4*x-15;
    // g=4*pow(x,5)+4*pow(x,3)-7*pow(x,2)-2*pow(x,4)+x-12;
    
    // f=9*pow(x,5)+2*pow(x,4)*y*z-189*pow(x,3)*pow(y,3)*z+117*pow(x,3)*y*pow(z,2)+3*pow(x,3)-42*pow(x,2)*pow(y,4)*pow(z,2)
    //                 +26*pow(x,2)*pow(y,2)*pow(z,3)+18*pow(x,2)-63*x*pow(y,3)*z+39*x*y*pow(z,2)+4*x*y*z+6;
    // g=6*pow(x,6)-126*pow(x,4)*pow(y,3)*z+78*pow(x,4)*y*pow(z,2)+pow(x,4)*y+pow(x,4)*z+13*pow(x,3)
    //     -21*pow(x,2)*pow(y,4)*z-21*pow(x,2)*pow(y,3)*pow(z,2)+13*pow(x,2)*pow(y,2)*pow(z,2)+13*pow(x,2)*y*pow(z,3)
    //     -21*x*pow(y,3)*z+13*x*y*pow(z,2)+2*x*y+2*x*z+2;
    // g=-3*pow(y,8)*pow(d,2)-2*pow(y,5)*pow(d,5)-2*pow(d,10)+5*pow(y,8)*d+3*pow(y,5)*pow(d,4)
    //     -7*pow(y,4)*pow(d,5)+5*pow(y,3)*pow(d,6)+10*pow(d,6)+pow(d,5)+3*pow(y,3)+4;
    // f=-9*pow(x,8)*pow(z,2)+6*pow(x,5)*pow(z,5)+2*pow(x,4)*pow(z,6)+9*pow(x,2)*z-9;
    // g=-7*pow(x,8)*pow(z,2)+2*pow(x,6)*y*pow(z,3)-3*pow(x,5)*pow(y,3)*pow(z,2)+pow(x,5)*y*pow(z,4)+10*pow(x,5)*pow(z,5)-9*pow(x,4)*pow(y,3)*pow(z,3)-pow(x,2)*pow(y,5)*pow(z,3)-6*pow(x,2)*y*pow(z,7)
    //   -2*x*pow(z,9)+2*pow(x,5)*pow(y,4)-8*pow(x,2)*pow(y,3)*pow(z,4)-pow(x,2)*pow(z,7)-4*x*pow(y,5)*pow(z,2)+3*x*pow(y,4)*pow(z,2)-2*x*pow(y,2)*pow(z,4)-7*pow(x,5)*z+7*pow(x,3)*pow(y,3)-2*x*pow(y,2)*pow(z,3)+6*y*pow(z,5)+9*pow(x,2)*pow(y,3)-8*x*pow(z,4)+9*pow(x,3)*z+2*pow(x,2)*pow(z,2)+5;
    // f=-7*pow(x,6)*y*pow(z,3)-8*pow(x,4)*pow(y,5)*z+10*pow(x,3)*pow(y,5)*pow(z,2)-3*pow(x,2)*pow(y,3)*pow(z,5)+pow(x,2)*pow(y,2)*pow(z,6)-2*x*y*pow(z,8)+10*pow(y,5)*pow(z,5)+pow(x,6)*y*pow(z,2)
    //    +8*pow(x,5)*pow(y,4)-9*pow(y,7)*pow(z,2)+pow(y,5)*pow(z,4)-pow(x,5)*y*pow(z,2)+10*pow(x,4)*pow(z,4)-5*pow(x,2)*pow(y,6)-3*x*pow(y,5)*pow(z,2)-6*pow(x,2)*pow(y,3)*pow(z,2)-x*pow(y,6)-5*pow(y,5)*pow(z,2)
    //   -9*pow(z,7)+2*pow(x,4)*pow(y,2)+6*pow(x,4)*y*z+7*x*pow(z,5)-pow(y,4)*pow(z,2)+9*pow(x,4)*z+4*pow(x,2)*pow(y,3)-9*x*pow(y,3)*z-3*x*y*pow(z,3)+9*y*pow(z,4)-10*pow(x,4)-6*x*pow(z,3)-9*y*pow(z,3)+6*x*y*z+8;
    // g=6*pow(x,10)+9*pow(x,8)*pow(y,2)-2*pow(y,4)*pow(z,6)-9*pow(x,4)*pow(z,5)+2*pow(x,3)*pow(y,4)*pow(z,2)+5*pow(x,8)+9*pow(x,6)*pow(z,2)+pow(x,3)*pow(y,5)-2*pow(x,3)*pow(y,2)*pow(z,3)+6*pow(x,2)*pow(z,6)+4*pow(y,6)*pow(z,2)+5*pow(y,4)*pow(z,4)-pow(x,7)-10*pow(x,6)*z+10*pow(x,3)*pow(y,2)*pow(z,2)+8*pow(y,6)*z-4*pow(y,5)*pow(z,2)+10*pow(y,2)*pow(z,5)-5*pow(x,2)*pow(y,2)*pow(z,2)+2*pow(x,2)*y*pow(z,3)-7*x*pow(y,3)*pow(z,2)-10*pow(x,4)+5*pow(y,4)+8;
    // f=-2*pow(x,5)*pow(y,2)*pow(z,3)-3*x*pow(y,7)*pow(z,2)+10*x*pow(y,5)*pow(z,4)-4*pow(y,10)+4*pow(y,6)*pow(z,4)+7*pow(x,8)*z-5*pow(x,5)*pow(y,4)+4*pow(x,4)*pow(z,5)+7*x*pow(y,2)*pow(z,6)-4*pow(y,5)*pow(z,4)+6*pow(x,6)*y*z+7*pow(x,5)*pow(y,2)*z+5*pow(x,5)*pow(z,3)-7*pow(x,4)*pow(y,3)*z-10*pow(x,2)*pow(z,6)+2*x*pow(y,5)*pow(z,2)-2*pow(x,5)*y*z-10*pow(x,3)*pow(y,3)*z+x*pow(y,4)*z-5*pow(y,4)*pow(z,2)+9*pow(y,3)*pow(z,2)+2*pow(y,2)*pow(z,3)+8*pow(x,3)*y-5*x*y*pow(z,2)+10;
    // // clpoly::polynomial_<clpoly::ZZ,clpoly::lex> f_,g_;
    // // clpoly::poly_convert(f,f_);
    // // clpoly::poly_convert(g,g_);
    // std::cout<<"f:"<<f<<std::endl;
    // std::cout<<"g:"<<g<<std::endl;
    // // std::cout<<clpoly::polynomial_GCD(f*g,g*g)<<std::endl;
    // clpoly::lex_<clpoly::custom_var_order> mo(std::vector<clpoly::variable>({z,y,x}));
    // clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> p1(&mo);
    // clpoly::poly_convert(g,p1);
    // std::cout<< p1<<std::endl;
    clpoly::polynomial_QQ pq_1(1);
    clpoly::polynomial_QQ pq_2(2);
    std::cout<<pq_1<<" "<<pq_2<<" "<<pq_1/pq_2<<std::endl;
    return 0;
  

}