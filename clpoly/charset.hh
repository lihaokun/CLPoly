/**
 * @file charset.hh
 * @author ntimesp(nxp@mail.ustc.edu.cn) 李昊坤(ker@pm.me)  
 * @brief  charset
 * 
 * 
 */

#ifndef CLPOLY_CHARSET_HH
#define CLPOLY_CHARSET_HH
#include <clpoly/polynomial.hh>

namespace clpoly{
    template <class Tc>
    bool rankcompare
    (const polynomial_<Tc,lex_<custom_var_order>> &F,
    const polynomial_<Tc,lex_<custom_var_order>> &G)
    {
        assert(comp_consistent(F.comp(),G.comp()));
        //比较两个多项式秩的大小 rank(F)<rank(G)?
        if(is_number(G))    return false;
        else if(is_number(F))   return true;
        //都不是常数多项式
        // auto m=F.comp().comp.v_map;
        // auto x=get_last_var_deg(F);
        // auto y=get_last_var_deg(G);
        auto & comp=F.comp();
        auto x=get_first_var(F);
        auto y=get_first_var(G);
        
        if(comp(y,x)) return true;
        if(x!=y) return false;
        if(get_first_deg(F)<get_first_deg(G))    return true;
        return false;
    }

    //判断F是否对G约化
    template <class Tc>
    bool is_reduced
    (const polynomial_<Tc,lex_<custom_var_order>> &F,
    const polynomial_<Tc,lex_<custom_var_order>> &G)
    {
        if(F.empty()) return true;
        //判断F是否对G约化
        if(is_number(G)) return false;
        if(is_number(F)) return false;

        if(degree(F,get_first_var(G))<get_first_deg(G))  return true;
        else return false;
    }

    template <class Tc,class InputIt>
    std::vector<polynomial_<Tc,lex_<custom_var_order>>> basset
    (InputIt first,InputIt second)
    {
        if(first==second) return {};
        std::vector<polynomial_<Tc,lex_<custom_var_order>>> F(first,second),B,Temp;
        while(!F.empty())
        {
            //找最小秩的多项式b
            std::size_t min=0;
            for(std::size_t i=1;i<F.size();++i)
                if(rankcompare(F[i],F[min]))    min=i;
            B.push_back(F[min]);
            //求出新的F
            Temp.clear();
            for(auto i:F)
                if(is_reduced(i,F[min]))  Temp.push_back(i);
            F=std::move(Temp);
            
        }
        return B;
    }

    template <class Tc>
    std::vector<polynomial_<Tc,lex_<custom_var_order>>> charset
    (const std::vector<polynomial_<Tc,lex_<custom_var_order>>>& P)
    {
        if(P.size()==0) return P;
        std::vector<polynomial_<Tc,lex_<custom_var_order>>> R,C;
        std::set<polynomial_<Tc,lex_<custom_var_order>>> F(P.begin(),P.end());
        while(true)
        {
            R.clear();
            C=basset<Tc>(F.begin(),F.end());
            //std::cout<<"C="<<C<<std::endl;
            // SHOW(C)
            if(is_number(C[0])) return C;
            for(auto f:F)
            {
                for(auto c=C.rbegin();c<C.rend();c++)
                {
                    f=prem(f,*c,get_first_var(*c),false);
                }
                if(!f.empty()) 
                {   
                    R.push_back(f);
                
                }
            }

            //std::cout<<"R="<<R<<std::endl;
            // SHOW(R)
            if (R.empty()) break;
            for (auto &i:R)
                F.insert(std::move(i));
            //std::cout<<"F="<<F<<std::endl;
            
        }
        return C;
    }
}

#endif