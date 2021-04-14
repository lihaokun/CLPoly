#include <clpoly/clpoly.hh>

int main()
{
    clpoly::variable x1("x1"),x2("x2"),x3("x3"),x4("x4");
    clpoly::lex_<clpoly::custom_var_order> mo(clpoly::custom_var_order({x4,x3,x2,x1}));
    clpoly::polynomial_ZZ T1,T2,T3,T4,f;
    T1=(x1+1)*pow(x1,2);//pow(x1,4)-pow(x1,3);
    T2=x2;//pow(x2,2)+(1+x1)*x2;
    T3=(x3-5)*(x3-4)*(x2-3)*(x2-2)*(x1-1);
    f=x1+x2;//x1*(x1-x2)*(x1+x2);

    clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>> T1_(&mo),T2_(&mo),T3_(&mo),f_(&mo);
    std::vector<clpoly::polynomial_<clpoly::ZZ,clpoly::lex_<clpoly::custom_var_order>>> T;
    clpoly::poly_convert(T1,T1_);
    clpoly::poly_convert(T2,T2_);
    clpoly::poly_convert(T3,T3_);
    clpoly::poly_convert(f,f_);
    T.push_back(T1_);
    T.push_back(T2_);
    //T.push_back(T3_);

    //std::cout<<"T="<<T<<std::endl;
    auto W=clpoly::wrsd(T,f_);
    std::cout<<"H="<<std::endl;
    for(auto i:W.first)
        std::cout<<i<<std::endl;
    std::cout<<"G="<<std::endl;
    for(auto i:W.second)
        std::cout<<i<<std::endl;

    return 0;
}