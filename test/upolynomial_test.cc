#include<clpoly/upolynomial.hh>

int main(int argc, char const *argv[])
{
    clpoly::upolynomial_<clpoly::ZZ> p={{4,2},{2,1},{0,2}};
    clpoly::upolynomial_<clpoly::ZZ> p2={{5,2},{2,-2},{1,1}};
    std::cout<<p*p2<<std::endl;
    return 0;
}
