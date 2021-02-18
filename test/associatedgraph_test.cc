#include <clpoly/clpoly.hh>
#include <iostream>

int main(int argc, char const *argv[])
{
    
    int n=5;
    std::vector<clpoly::polynomial_ZZ> F;
    for (int i=0;i<n;++i)
    {
        clpoly::variable xk("x"+std::to_string(i+1));
        clpoly::variable xk1("x"+std::to_string(i+2));
        clpoly::variable xk2("x"+std::to_string(i+3));
        clpoly::variable xk3("x"+std::to_string(i+4));
        
        clpoly::polynomial_ZZ f1=xk*xk3-xk1*xk2;
        F.push_back(f1);
    }
    auto G3=clpoly::associatedgraph(F);
    std::cout<<clpoly::graph_diff_score(chordal_completion(G3),G3)<<std::endl;
    std::cout<<clpoly::graph_diff_score(connected_branch_graph(G3),G3)<<std::endl; 
    std::vector<clpoly::variable> l={
       clpoly::variable("x1"),
       clpoly::variable("x2"),
       clpoly::variable("x3"),
       clpoly::variable("x8"),
       clpoly::variable("x7"),
       clpoly::variable("x6"),
       clpoly::variable("x4"),
       clpoly::variable("x5")
    };

    std::cout<<clpoly::elimination_height(G3,l)<<std::endl;
    auto md=__polynomial_m_d(F,l);
    std::cout<<md.first<<" "<<md.second<<std::endl;

    
    l=std::vector<clpoly::variable>();
    for (auto index=1;index<=10;++index)
        l.push_back(clpoly::variable("x"+std::to_string(index)));
    std::cout<<l<<std::endl;
    F=clpoly::random_polynomials<clpoly::ZZ>(l,2,0.2,0.5,10,-10);
    std::cout<<F<<std::endl;
    auto G=clpoly::associatedgraph(F);
    std::cout<<chordal_completion(G)<<std::endl;
    l=clpoly::chordal_completion(G3,G);
    std::cout<<l<<std::endl;
    std::cout<<G3<<std::endl;
    std::cout<<clpoly::graph_diff_score(G3,G)<<std::endl;
    std::cout<<clpoly::graph_diff_score(clpoly::connected_branch_graph(G),G)<<std::endl;
    return 0;
}
