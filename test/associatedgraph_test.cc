#include <clpoly/clpoly.hh>
#include <iostream>

int main(int argc, char const *argv[])
{
    
    // int n=5;
    // std::vector<clpoly::polynomial_ZZ> F;
    // for (int i=0;i<n;++i)
    // {
    //     clpoly::variable xk("x"+std::to_string(i+1));
    //     clpoly::variable xk1("x"+std::to_string(i+2));
    //     clpoly::variable xk2("x"+std::to_string(i+3));
    //     clpoly::variable xk3("x"+std::to_string(i+4));
        
    //     clpoly::polynomial_ZZ f1=xk*xk3-xk1*xk2;
    //     F.push_back(f1);
    // }
    // auto G3=clpoly::associatedgraph(F);
    // std::cout<<clpoly::graph_diff_score(chordal_completion(G3),G3)<<std::endl;
    // std::cout<<clpoly::graph_diff_score(connected_branch_graph(G3),G3)<<std::endl; 
    // std::vector<clpoly::variable> l={
    //    clpoly::variable("x1"),
    //    clpoly::variable("x2"),
    //    clpoly::variable("x3"),
    //    clpoly::variable("x8"),
    //    clpoly::variable("x7"),
    //    clpoly::variable("x6"),
    //    clpoly::variable("x4"),
    //    clpoly::variable("x5")
    // };

    // std::cout<<clpoly::elimination_height(G3,l)<<std::endl;
    // auto md=__polynomial_m_d(F,l);
    // std::cout<<md.first<<" "<<md.second<<std::endl;

    
    // l=std::vector<clpoly::variable>();
    // for (auto index=1;index<=10;++index)
    //     l.push_back(clpoly::variable("x"+std::to_string(index)));
    // std::cout<<l<<std::endl;
    // F=clpoly::random_polynomials<clpoly::ZZ>(l,2,0.2,0.5,10,-10);
    // std::cout<<F<<std::endl;
    // auto G=clpoly::associatedgraph(F);
    // std::cout<<chordal_completion(G)<<std::endl;
    // l=clpoly::chordal_completion(G3,G);
    // std::cout<<l<<std::endl;
    // std::cout<<G3<<std::endl;
    // std::cout<<clpoly::graph_diff_score(G3,G)<<std::endl;
    // std::cout<<clpoly::graph_diff_score(clpoly::connected_branch_graph(G),G)<<std::endl;
    
    // auto a=clpoly::variable("a");
    // auto t=clpoly::variable("t");
    // auto x=clpoly::variable("x");
    // auto y=clpoly::variable("y");
    // auto h=clpoly::variable("h");
    // auto s=clpoly::variable("s");
    // auto b=clpoly::variable("b");
    // auto c=clpoly::variable("c");
    // auto r=clpoly::variable("r");
    
    // // clpoly::polynomial_ZZ f1=t*t-2*x*t+x*x+4*y*y-80*y+396;
    // // clpoly::polynomial_ZZ f2=5*a*a*t*t-2*a*t*x-8*a*t*y+x*x+4*y*y-4;
    // clpoly::polynomial_ZZ f0=a*a*h*h-4*s*(s-a)*(s-b)*(s-c);
    // clpoly::polynomial_ZZ f1=2*r*h-b*c;
    // clpoly::polynomial_ZZ f2=2*s-a-b-c;
    
    // clpoly::polynomial_ZZ f3=a;
    // clpoly::polynomial_ZZ f4=b;
    // clpoly::polynomial_ZZ f5=c;
    // clpoly::polynomial_ZZ f6=r;
    // clpoly::polynomial_ZZ f7=h;
    // clpoly::polynomial_ZZ f10=a+b-c;
    // clpoly::polynomial_ZZ f11=b+c-a;
    // clpoly::polynomial_ZZ f12=-b+c+a;
    
    
    
    // std::vector<clpoly::polynomial_ZZ> F={
    //     f0,
    //     f1,
    //     f2,f3,f4,f5,f6,f7,f10,f11
    // };
    // std::cout<<F;
    
    // auto G3=clpoly::associatedgraph(F);
    // std::cout<<G3;
    
    int n=9;
    auto l=std::vector<clpoly::variable>();
    for (auto index=1;index<=n;++index)
        l.push_back(clpoly::variable("x["+std::to_string(index)+"]"));
    std::cout<<l<<std::endl;
    auto F=std::vector<clpoly::polynomial_ZZ>();
    for(auto index=0;index<n-3;++index)
        F.push_back(l[index]*l[index+3]-l[index+1]*l[index+2]);
    std::cout<<F<<std::endl;
    auto p=clpoly::peo(F);
    std::cout<<p<<std::endl;
    return 0;
}
