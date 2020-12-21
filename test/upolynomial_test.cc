#include<clpoly/upolynomial.hh>
#include <clpoly/realroot.hh>
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
    return 0;
}
