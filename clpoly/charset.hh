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

    template <class Tc>
    std::pair<
        std::vector<std::vector<polynomial_<Tc,lex_<custom_var_order>>>>,
        std::vector<std::vector<polynomial_<Tc,lex_<custom_var_order>>>>
    > wrsd
    (
        const std::vector<polynomial_<Tc,lex_<custom_var_order>>> &T,
        const polynomial_<Tc,lex_<custom_var_order>> &f
    )
    {
        std::vector<std::vector<polynomial_<Tc,lex_<custom_var_order>>>> H,G;

        // std::cout<<"T="<<T<<std::endl;
        // std::cout<<"f="<<f<<std::endl;

        auto g=f;
        for(auto c=T.rbegin();c<T.rend();c++)
            if(!is_reduced(g,*c))
                g=prem(g,*c,get_first_var(*c));  
        if(f!=g) return wrsd(T,g);

        if(f.empty()) 
        {
            H.push_back(T);
            return std::make_pair(H,G);
        }
        if(is_number(f)) 
        {
            G.push_back(T);
            return std::make_pair(H,G);
        }
        //最后出现变量为x_n(数组下标n-1)
        auto m=T[0].comp().comp.v_map;
        auto n=m.size()+1-m[get_first_var(f)];
        auto x=get_first_var(f);
        //std::cout<<"n="<<n<<std::endl;

        auto res=f;
        for(int i=n-1;i>=0;i--)
        {   
            //std::cout<<"res="<<res<<std::endl;
            res=resultant(res,T[i],get_first_var(T[i]));
        }
        //std::cout<<"res="<<res<<std::endl;
        if(!res.empty())
        {
            G.push_back(T);
            return std::make_pair(H,G);
        }

        //计算正规子结式链
        auto subreschain=subresultant(T[n-1],f,x);
        std::vector<polynomial_<Tc,lex_<custom_var_order>>> regularsubreschain;
        for(auto i=0;i<subreschain.size();i++)
            //由于deg(f)<deg(T[n-1])（因为f关于T约化），所以T[n-1]自动被这个循环加进去了，不用补在最后
            if(degree(subreschain[i],x)==i)
                regularsubreschain.push_back(subreschain[i]);
        //auto v=regularsubreschain.size()-2;

        // std::cout<<"subreschain="<<subreschain<<std::endl;
        // std::cout<<"regularsubreschain="<<regularsubreschain<<std::endl;

        if(n==1)
        {
            auto q=pquo(T[n-1],regularsubreschain[1],x);
            std::vector<polynomial_<Tc,lex_<custom_var_order>>> TT;
            TT.push_back(q);

            //std::cout<<"q="<<q<<std::endl;

            auto W=wrsd(TT,f);
            TT[0]=regularsubreschain[1];
            H.push_back(TT);
            G=W.second;
        }
        else
        {
            std::vector<polynomial_<Tc,lex_<custom_var_order>>> TT(T.begin(),T.begin()+n-1);
            auto W=wrsd(TT,regularsubreschain[0]);
            auto Hi=W.first;
            auto Gi=W.second;
            for(auto t:Gi) 
            {
                t.push_back(T[n-1]);
                G.push_back(t);
            }
            auto i=0;
            while(!Hi.empty())
            {
                i++;
                std::vector<std::vector<polynomial_<Tc,lex_<custom_var_order>>>> Hii,Gii;
                auto R=leadcoeff(regularsubreschain[i]);

                // std::cout<<"S_d="<<regularsubreschain[i]<<std::endl;
                // std::cout<<"R="<<R<<std::endl;

                for(auto h:Hi)
                {
                    W=wrsd(h,R);
                    Hii.insert(Hii.end(),W.first.begin(),W.first.end());
                    Gii.insert(Gii.end(),W.second.begin(),W.second.end());
                }
                for(auto g:Gii)
                {
                    g.push_back(regularsubreschain[i]);
                    H.push_back(g);
                    auto q=pquo(T[n-1],regularsubreschain[i],x);
                    if(degree(q,x)>0)
                    {
                        g.erase(g.end());
                        g.push_back(q);
                        W=wrsd(g,f);
                        G.insert(Gii.end(),W.second.begin(),W.second.end());
                    }
                }
                Hi=Hii;
                Gi=Gii;
            }
        }

        //需要补上n开始的最后几个多项式
        for(auto &i:H)  i.insert(i.end(),T.begin()+n,T.end());
        for(auto &i:G)  i.insert(i.end(),T.begin()+n,T.end());

        // std::cout<<"H="<<std::endl;
        // for(auto i:H)
        //     std::cout<<i<<std::endl;
        // std::cout<<"G="<<std::endl;
        // for(auto i:G)
        //     std::cout<<i<<std::endl;

        return std::make_pair(H,G);
    }
}

#endif