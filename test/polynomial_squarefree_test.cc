#include <clpoly/clpoly.hh>

int main(int argc, char const *argv[])
{
    clpoly::variable x1("x1"),x2("x2"),x3("x3"),x4("x4");
    clpoly::lex_<clpoly::custom_var_order> mo(clpoly::custom_var_order({x4,x3,x2,x1}));
    clpoly::polynomial_ZZ f,g;
    f=-72*x2*pow(x1,8)+360*x2*pow(x1,6)+432*pow(x1,7)+80*pow(x1,6)-2160*pow(x1,5)-400*pow(x1,4);
    //g=(x2-2)*(x1-1);
    clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> f_(&mo),g_(&mo);
    clpoly::poly_convert(f,f_);
    //clpoly::poly_convert(g,g_);
    auto s=clpoly::squarefreefactorize(f_);
    std::cout<<f_<<std::endl;
    for(auto i:s)
    {
        std::cout<<"["<<i.first<<","<<i.second<<"]  ";
        // std::cout<<"["<<i<<"]  ";
    }
}