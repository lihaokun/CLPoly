#include<clpoly/upolynomial.hh>
#include <clpoly/realroot.hh>
#include <clpoly/clpoly.hh>
int main(int argc, char const *argv[])
{
    clpoly::upolynomial_<clpoly::ZZ> p={{4,2},{2,1},{0,2}};
    clpoly::upolynomial_<clpoly::ZZ> p2={{5,2},{2,-2},{1,1}};
    std::cout<<p*p2<<std::endl;

    clpoly::upolynomial_<clpoly::ZZ> f={{0,1}};
    for (int i=1;i<21;++i)
    {
        p={{1,1},{0,i}};
        f=f*p;
    }
    p={{19,-1}};
    p2={{0,clpoly::pow(clpoly::ZZ(10),9)}};

    f=f*p2+p;
    
    std::cout<<f<<std::endl;
    auto l=clpoly::uspensky(f);
    for (auto &i:l)
    {
        std::cout<<"{"<<i.first<<","<<i.second<<"}"<<" ";
    }
    std::cout<<std::endl;
    auto x=clpoly::variable("x");
    clpoly::polynomial_ZZ Fx;
    for(auto &i:f)
    {
        Fx.push_back({clpoly::monomial({{x,i.first.deg()}}),i.second});
    }
    // Fx.normalization();
    std::cout<<Fx<<std::endl;
    auto roots=clpoly::realroot<clpoly::grlex>({Fx,Fx*Fx,x*x-2,x*x});
    for (auto &i:roots.first)
    {
        std::cout<<"{"<<i.first<<","<<i.second<<"}"<<" ";
    }
    std::cout<<std::endl;
    for(auto &i:roots.second)
    {
        std::cout<<"{";
        for(auto &j:i)
            std::cout<<j.first<<"*"<<j.second<<" ";
        std::cout<<"} ";
    }
    std::cout<<std::endl;
    // auto Rp1=clpoly::random_polynomial<clpoly::ZZ>({x},100,0.1,10,-10);
    //  clpoly::upolynomial_<clpoly::ZZ> f1(Rp1);
    //  f1=Rp1;
    // auto x=clpoly::variable("x");
    // for (auto  i=1;i<=100;++i)
    // {
    //     auto Rp1=clpoly::random_polynomial<clpoly::ZZ>({x},100,0.1,10,-10);
    //     auto Rp2=clpoly::random_polynomial<clpoly::ZZ>({x},100,0.1,10,-10);
    //     // std::cout<<Rp1<<std::endl;
    //     // std::cout<<Rp2<<std::endl;
        
    //     clpoly::upolynomial_<clpoly::ZZ> f1;
    //     clpoly::upolynomial_<clpoly::ZZ> g1;
    //     clpoly::poly_convert(Rp1,f1);
    //     clpoly::poly_convert(Rp2,g1);
        
    //     std::cout<<"f["<<i<<"]="<<f1<<std::endl;
    //     std::cout<<"g["<<i<<"]="<<g1<<std::endl;
    //     std::cout<<"o["<<i<<"]="<<clpoly::polynomial_GCD(f1*f1,g1*f1)<<std::endl;
    // }
    return 0;
}
