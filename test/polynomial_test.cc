#include <clpoly.hh>
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
    clpoly::variable r("r");
    
    // clpoly::polynomial_ZZ p=2*y*z*x*x*z+d*d*d*y*z*z+1;
    // std::cout<< p<<std::endl;
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
    // std::cout<<"4*x^4*y^2*z^4+4*x^2*y^2*d^3*z^4+y^2*d^6*z^4+4*x^2*y*z^2+2*y*d^3*z^2+1\n";
    // std::cout<< p<<std::endl;
    // clpoly::polynomial_ZZ f=x*pow(y,2)+1;
    // clpoly::polynomial_ZZ g=2*pow(y,3)-pow(y,2)+pow(x,2)*y;

    clpoly::polynomial_ZZ f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},10,0.2,10,-10);
    clpoly::polynomial_ZZ g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z,d},10,0.2,10,-10);
    clpoly::polynomial_ZZ o;
     clpoly::polynomial_ZZ o1;
    // std::cout<<"g:="<<g<<":"<<std::endl;
    // std::cout<<"f:="<<f<<":"<<std::endl;
    auto t=clock();
    double t1;
    // std::cout<<"o:="<<clpoly::prem(g,f,y)<<":"<<std::endl;
    // std::cout<<"t:="<<(double(clock()-t)/CLOCKS_PER_SEC)<<";"<<std::endl;
    // std::cout<<"st := time():o1:= expand(prem(g, f, y)):time() - st;o1-o;"<<std::endl;
    
    // f=-pow(x,2)*pow(z,3) - pow(x,4) - pow(z,4) + pow(x,2) + 2*pow(z,2) - 1;
    // g=-pow(r,2)*pow(x,2) + pow(x,4) + pow(x,2)*pow(z,2) + pow(z,4) - 2*pow(z,2) + 1;
    //g=-4*pow(x,10)+8*pow(x,9)*z-7*pow(x,8)*y*z-4*pow(x,7)*pow(z,3)+10*pow(x,4)*pow(y,2)*pow(z,4)-pow(x,4)*pow(z,6)-7*pow(x,3)*pow(y,6)*z+2*pow(x,3)*y*pow(z,6)-2*x*pow(y,4)*pow(z,5)-4*pow(y,8)*pow(z,2)-6*pow(y,7)*pow(z,3)+2*y*pow(z,9)-6*pow(x,7)*y*z+7*pow(x,5)*pow(y,3)*z-2*pow(x,3)*pow(y,2)*pow(z,4)-7*pow(x,2)*pow(y,5)*pow(z,2)-pow(x,2)*pow(y,2)*pow(z,5)-8*x*pow(y,4)*pow(z,4)-8*x*y*pow(z,7)-10*pow(y,8)*z+8*pow(x,8)+5*pow(x,7)*y+6*pow(x,7)*z-5*pow(x,5)*pow(y,3)-5*pow(x,4)*pow(y,4)+9*pow(x,3)*pow(y,3)*pow(z,2)-9*x*pow(y,6)*z+x*pow(y,4)*pow(z,3)-7*y*pow(z,7)-10*pow(x,5)*pow(z,2)+9*pow(x,3)*pow(y,4)+3*pow(x,3)*pow(y,2)*pow(z,2)+6*pow(x,2)*pow(y,3)*pow(z,2)-8*x*pow(y,4)*pow(z,2)+7*pow(y,2)*pow(z,5)-4*y*pow(z,6)+8*pow(x,4)*pow(y,2)-pow(x,4)*y*z-pow(x,2)*y*pow(z,3)+8*pow(x,2)*pow(z,4)-2*x*pow(y,5)-8*pow(y,6)+2*y*pow(z,5)+pow(x,2)*pow(y,2)*z+pow(x,2)*y*pow(z,2)+4*x*y*pow(z,3)+3*pow(y,4)*z-9*pow(x,4)+7*pow(x,3)*z+3*x*pow(z,3)+4*pow(y,2)*pow(z,2)-4*pow(z,4)-6*pow(x,2)*y+8*x*pow(y,2)+7*x*y*z-8*pow(y,3)-2*x*y-3*x*z-8*pow(z,2)+1;
    //f=-8*pow(x,9)*z-pow(x,8)*pow(y,2)-2*pow(x,6)*pow(y,4)-6*pow(x,5)*pow(z,5)+2*pow(x,3)*pow(y,6)*z-3*pow(x,3)*pow(y,3)*pow(z,4)+2*pow(x,2)*pow(y,4)*pow(z,4)+9*pow(x,2)*pow(z,8)+3*x*pow(y,9)+3*x*pow(y,7)*pow(z,2)+5*pow(y,10)-5*pow(x,9)+5*pow(x,7)*pow(z,2)-6*pow(x,5)*pow(z,4)-3*pow(x,3)*pow(y,6)+pow(x,3)*pow(y,4)*pow(z,2)-2*pow(x,2)*pow(y,3)*pow(z,4)-4*pow(x,2)*pow(y,2)*pow(z,5)-8*pow(x,2)*y*pow(z,6)-3*x*pow(y,6)*pow(z,2)-9*x*pow(y,3)*pow(z,5)+2*pow(x,3)*pow(y,5)-pow(x,2)*pow(y,6)+4*x*pow(y,5)*pow(z,2)+2*pow(y,4)*pow(z,4)-3*pow(x,2)*pow(y,5)+pow(x,2)*pow(y,4)*z-8*pow(x,2)*pow(y,2)*pow(z,3)+8*pow(x,2)*pow(z,5)+5*x*pow(z,6)+7*pow(x,4)*pow(z,2)+5*pow(x,3)*y*pow(z,2)-7*pow(x,2)*pow(y,2)*pow(z,2)-5*x*pow(y,3)*pow(z,2)-4*x*y*pow(z,4)+3*pow(x,4)*y-3*pow(x,4)*z+pow(x,3)*y*z+6*pow(x,3)*pow(z,2)+9*pow(x,2)*pow(y,2)*z-2*x*pow(y,4)-9*x*pow(z,4)-6*pow(x,2)*y*z+5*pow(x,2)*pow(z,2)+9*x*y*pow(z,2)-8*pow(y,3)*z-pow(y,2)*pow(z,2)-9*x*pow(z,2)+y*z+z;
    //g=-6*pow(x,10)+5*pow(x,9)*y+2*pow(x,4)*pow(y,2)*pow(z,4)+10*pow(x,3)*pow(y,6)*z+2*pow(x,3)*pow(z,7)-5*pow(x,2)*pow(y,5)*pow(z,3)-7*x*pow(y,6)*pow(z,3)+7*pow(y,10)+9*pow(y,8)*pow(z,2)-3*pow(y,6)*pow(z,4)+8*pow(y,2)*pow(z,8)+3*pow(x,8)*y-2*pow(x,3)*pow(y,6)+5*pow(x,3)*pow(y,4)*pow(z,2)-3*pow(x,2)*pow(y,7)+3*x*pow(y,8)-2*x*pow(y,3)*pow(z,5)+8*pow(y,8)*z+pow(x,7)*y+3*pow(x,4)*pow(y,4)+9*pow(x,4)*pow(y,2)*pow(z,2)+6*pow(x,3)*pow(y,4)*z-7*pow(x,2)*pow(y,2)*pow(z,4)+5*x*pow(y,7)+4*x*pow(y,3)*pow(z,4)+9*pow(y,7)*z+3*pow(y,6)*pow(z,2)-4*pow(y,3)*pow(z,5)-4*pow(x,5)*pow(y,2)-6*pow(x,3)*pow(y,3)*z+6*pow(x,2)*pow(y,3)*pow(z,2)-pow(x,2)*pow(y,2)*pow(z,3)+7*x*pow(y,4)*pow(z,2)-9*x*y*pow(z,5)+7*pow(y,5)*pow(z,2)+pow(x,2)*pow(z,4)-6*pow(x,3)*pow(y,2)+2*pow(x,2)*y*pow(z,2)-10*x*pow(y,3)*z-5*pow(y,2)*pow(z,3)-6*pow(z,5)-10*pow(x,2)*pow(z,2)-4*x*y*pow(z,2)-3*x*pow(z,3)-10*pow(z,4)-7*pow(x,3)-5*pow(x,2)*y-3*pow(y,3)+y*pow(z,2)+4*pow(z,3)-7*x*y+x*z-9;
    //f=3*pow(x,9)*z+4*pow(x,7)*pow(y,2)*z+6*pow(x,6)*pow(y,4)-8*pow(x,4)*pow(y,4)*pow(z,2)-8*pow(x,3)*y*pow(z,6)+pow(x,2)*pow(y,6)*pow(z,2)+9*x*pow(y,9)-x*pow(y,8)*z+4*x*pow(y,5)*pow(z,4)+6*pow(x,3)*pow(y,4)*pow(z,2)-9*pow(x,3)*pow(y,2)*pow(z,4)-3*pow(x,3)*y*pow(z,5)-7*pow(x,3)*pow(z,6)-4*pow(x,2)*pow(y,7)+9*pow(x,2)*pow(y,6)*z-4*x*pow(y,7)*z-2*pow(y,6)*pow(z,3)+6*pow(y,3)*pow(z,6)-5*pow(y,2)*pow(z,7)+7*pow(x,4)*y*pow(z,3)-pow(x,3)*y*pow(z,4)-4*pow(x,2)*pow(y,2)*pow(z,4)+4*x*pow(y,6)*z+x*pow(y,5)*pow(z,2)-10*x*y*pow(z,6)-3*pow(y,7)*z+2*pow(x,6)*y+2*pow(x,4)*y*pow(z,2)-7*pow(x,2)*y*pow(z,4)+3*x*pow(y,2)*pow(z,4)+2*x*pow(z,6)+8*pow(y,6)*z+pow(y,5)*pow(z,2)-4*pow(x,4)*pow(z,2)+pow(y,5)*z+5*pow(y,3)*pow(z,3)-10*pow(y,2)*pow(z,4)-8*y*pow(z,5)-4*pow(x,3)*y*z+7*x*pow(y,4)-2*pow(y,4)*z-3*pow(y,3)*pow(z,2)-pow(y,2)*pow(z,3)-7*x*y*pow(z,2)+10*pow(y,3)*z-5*pow(z,4)+6*x*pow(y,2)+6*pow(z,3)-3*pow(z,2)-10*y+2*z-7;
    // g=4*x*pow(y,3)*z-4*x*pow(y,2)*pow(z,2)-8*x*pow(z,4)-8*pow(y,4)*z-3*pow(x,2)*y*z+2*pow(x,2)*pow(z,2)-8*x*pow(y,3)-x*pow(z,3)+2*pow(x,2)*z-5*x*pow(y,2)-8*x*pow(z,2)+9*y*pow(z,2)-4*pow(z,3)+5*x*y+3;
    // f=7*pow(x,2)*pow(y,3)+x*pow(y,4)-8*x*pow(z,4)+9*pow(y,5)-9*y*pow(z,4)+2*pow(x,4)-7*pow(x,2)*pow(y,2)-7*x*pow(y,3)+6*pow(y,4)-9*pow(y,2)*pow(z,2)-3*y*pow(z,3)+2*pow(z,4)+pow(x,2)*z+5*x*pow(y,2)-3*x*pow(z,2)-6*pow(y,3)+9*y*z+4;
    // f=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.2,10,-10);
    // g=clpoly::random_polynomial<clpoly::ZZ>({x,y,z},10,0.2,10,-10);
    //g=-6*pow(x,3)*pow(z,2)+4*pow(x,2)*pow(z,3)-7*pow(z,5)+pow(x,3)*y+10*x*pow(y,3)+3*pow(y,4)+y*pow(z,3)+2*pow(z,2)+x-3*z-3;
    //f=-6*pow(x,5)-10*pow(x,3)*pow(y,2)-9*pow(x,3)*y*z-8*pow(x,3)*pow(z,2)+4*pow(x,2)*pow(y,3)+5*pow(x,2)*y*pow(z,2)-7*x*y*pow(z,3)-3*pow(y,5)-9*pow(x,2)*pow(z,2)+9*y*pow(z,2)-3*pow(y,2)-7;
    // g=5*pow(x,5)+9*pow(x,4)*y+pow(x,4)*z-10*pow(x,3)*y*z+10*pow(x,3)*pow(z,2)-8*pow(x,2)*y*pow(z,2)-9*x*pow(z,4)+2*pow(y,3)*pow(z,2)-5*pow(x,4)+7*pow(x,3)*y-pow(x,2)*z+6*x*pow(z,2)+6*pow(y,3)-8;
    // f=-2*pow(x,3)*pow(z,2)+7*pow(x,2)*pow(y,2)*z+x*y*pow(z,3)+2*pow(y,5)-10*pow(z,5)-3*pow(x,3)*z-3*pow(y,3)*z+6*pow(x,3)-5*pow(y,3)-3*z+1;
    // g=5*pow(x,5)*pow(y,5)+7*pow(x,6)*pow(y,2)+8*pow(x,5)*pow(y,3)-6*pow(y,8)-2*pow(x,4)*y+3*x*pow(y,2)-6*pow(y,3)-3;
    // f=-4*pow(y,6)+7*pow(x,4)*y-10*pow(x,4)-4;
    // g=-5*pow(x,7)*y-3*pow(x,6)*y-4*pow(x,6)+9*x*pow(y,3);
    // f=3*pow(x,6)*pow(y,4)+pow(x,6)*pow(y,2)+3*pow(x,5)*pow(y,2)+7;

    // std::cout<<"t=Association[];t2=Association[];\n";
    // for (int i=0;i<1;++i)
    // {
    //     f=clpoly::random_polynomial<clpoly::ZZ>({x,y,d},10,0.1,10,-10);
    //     g=clpoly::random_polynomial<clpoly::ZZ>({x,y,d},10,0.1,10,-10);
    //     std::cout<<"g="<<g<<";"<<std::endl;
    //     std::cout<<"f="<<f<<";"<<std::endl;
    //     // f=-6*pow(x,4)-7*pow(x,2)*y-x*pow(y,3)*z-3*x*y-9;
    //     // g=10*pow(x,4)*y+2*pow(x,4)-7*pow(x,3)*pow(y,2)+9*pow(x,2)+9*pow(z,4)-9*y+2*z;
    //     // std::cout<<prem(f,g,x)<<std::endl;
    //     clpoly::__is_monomial_compression=true;
    //     t=clock();
    //     o=resultant(f,g,x);
    //     t1=(double(clock()-t)/CLOCKS_PER_SEC);
    //     std::cout<<"o="<<o<<";"<<std::endl;
    //     std::cout<<"t["<<i<<"]=AbsoluteTiming[o1=o-Expand[Resultant[f,g,x]]];"<<std::endl;
    //     printf("t2[%d]=%.6f;\n",i,t1);
    //     // std::cout<<"t2["<<i<<"]="<<std::fixed<<setprecision(6)<<t1<<";"<<std::endl;

    //     // clpoly::__is_monomial_compression=false;
    //     // t=clock();
    //     // o1=resultant(f,g,x);
    //     // t1=(double(clock()-t)/CLOCKS_PER_SEC);
    //     // std::cout<<"o="<<(o1==o)<<";"<<std::endl;
    //     // std::cout<<"time:"<<t1<<std::endl;
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

    // clpoly::Zp a(10,7);
    // clpoly::Zp b(-10,7);
    // a=-2;
    // std::cout<<a<<" "<<b<<" "<<a*b<<" "<<a/b<<" "<<-a<<std::endl;

    // g=-pow(r,2)*pow(x,2) + pow(x,4) + pow(x,2)*pow(z,2) + pow(z,4) - 2*pow(z,2) + 1;
    // std::cout<<g<<std::endl;
    // std::cout<<clpoly::polynomial_mod(std::move(g),7)<<std::endl;
    // a=std::move(b);
    // std::cout<<a<<b<<std::endl;
    f=pow(x,4)+25*pow(x,3)+145*pow(x,2)-171*x-360;
    g=pow(x,5)+14*pow(x,4)+15*pow(x,3)-pow(x,2)-14*x-15;
    // f=2*pow(x,4)-7*pow(x,3)-4*pow(x,2)-4*x-15;
    // g=4*pow(x,5)+4*pow(x,3)-7*pow(x,2)-2*pow(x,4)+x-12;
    std::cout<<clpoly::polynomial_GCD(f,g)<<std::endl;
    return 0;
}