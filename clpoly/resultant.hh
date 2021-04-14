/**
 * @file resultant.hh
 * @author 李昊坤(ker@pm.me)  ntimesp(nxp@mail.ustc.edu.cn)
 * @brief  定义结式运算 在变量小于4时用子子结式链 大于时用 Bezout 矩阵
 *         定义子结式运算
 * 
 * 
 */

#ifndef CLPOLY_RESULTANT_HH
#define CLPOLY_RESULTANT_HH

#include <clpoly/polynomial.hh>
#include <list>
#include <string>

namespace clpoly{
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> resultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        assert(comp_consistent(G.comp(),F.comp()));
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        resultant(O1,G1,F1,v);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(std::move(O1),O);
        return O;
    }
    
    
    template <class Tc>
    void resultant
        (   
            polynomial_<Tc,univariate_priority_order>&O,
            const polynomial_<Tc,univariate_priority_order>&F,
            const polynomial_<Tc,univariate_priority_order>&G,
            variable v
        )
    {
        O.clear();
        const univariate_priority_order &comp=G.comp();
       
        int64_t m=get_up_deg(F);
        int64_t l=get_up_deg(G);
        if (m<0 || l<0 || G.empty()||F.empty())
            return void();
       
        if (m==0)
        {
            O=pow(F,l);
            return void();
        }
        if (l==0)
        {
            O=pow(G,m);
            return  void();
        }
        
        if (m<l)
        {
            resultant(O,G,F,v);
            if ((l&1) &&(m&1))
                pair_vec_negate(O.data());
            return void();
        }
        std::list<std::pair<variable,int64_t>> vars;
        __pair_vec_variables(F.data(),vars);
        __pair_vec_variables(G.data(),vars);
        // std::cout<<vars.size()<<std::endl;
        if (vars.size()>3)
        {
            auto Bez =BezoutMatrix(F,G);
            O=det(Bez,m);
            return void();
        }

        polynomial_<Tc,univariate_priority_order> S_j_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_j(&comp);
        polynomial_<Tc,univariate_priority_order> S_r_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_r(&comp);
        
        polynomial_<Tc,univariate_priority_order> R_(&comp);
        polynomial_<Tc,univariate_priority_order> tmp1(&comp);
        polynomial_<Tc,univariate_priority_order> tmp2(&comp);
        int64_t j,r,u;
        // std::cout<<"F="<<F<<std::endl;
        // std::cout<<"G="<<G<<std::endl;
        
        if (l<m)
        {
            j=m-1;
            if(j==l)
            {
                prem(S_j,F,G,v);
                --j;
                //  std::cout<<j<<":"<<S_j<<std::endl;
                S_j_1=G;
            }
            else{
                leadcoeff(R_,G);
                if (j-l<2)
                    pair_vec_multiplies(S_r.data(),R_.data(),G.data(),comp);
                else 
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-l,comp);
                    pair_vec_multiplies(S_r.data(),tmp1.data(),G.data(),comp);
                }
                //  std::cout<<l<<":"<<S_r<<std::endl;
                prem(S_r_1,F,G,v);
                if (j-l & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<l-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=l-1;
            }

        }else
        {
            j=m;
            prem(S_j,F,G,v);
            pair_vec_negate(S_j.data());
            //  std::cout<<j-1<<":"<<S_j<<std::endl;
            if (!(--j))
            {
                O.data()=std::move(S_j.data());
                return void();
            }
            r=get_up_deg(S_j);
            if (r<0)
                return void();
            if (r<j)
            {
                
                if (r!=1)
                {
                    leadcoeff(R_,S_j);
                    if (j-r<2)
                        pair_vec_multiplies(S_r.data(),R_.data(),S_j.data(),comp);
                    else 
                    {
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_multiplies(S_r.data(),tmp1.data(),S_j.data(),comp);
                    }
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    O.data()=std::move(S_r.data());
                    return void();
                }
                prem(tmp1,G,S_j,v);
                //std::cout<<"G"<<":"<<G<<std::endl;
                //std::cout<<"S_"<<j<<":"<<S_j<<std::endl;
                
                //  std::cout<<r-1<<":"<<tmp1<<std::endl;
                leadcoeff(R_,G);
                pair_vec_div(S_r_1.data(),tmp1.data(),R_.data(),comp);
                if ((j-r) & 1==1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=r-1;
            }
            else
            {
                prem(tmp1,G,S_j,v);
                leadcoeff(R_,G);
                pair_vec_div(S_j_1.data(),tmp1.data(),R_.data(),comp);
                //  std::cout<<j-1<<":"<<S_j_1<<std::endl;
                --j;
                swap(S_j.data(),S_j_1.data());
            }
        }
        while (j>0)
        {
            r=get_up_deg(S_j);
            if (r<0)
                return void();
            if (r<j)
            {
                
                //std::cout<<"R_j"<<":"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                //std::cout<<"R_j+1"<<":"<<R_<<std::endl;
                if (r!=1)
                {
                    leadcoeff(tmp1,S_j);
                    if (j-r<2)
                    {
                        pair_vec_multiplies(tmp2.data(),tmp1.data(),S_j.data(),comp);
                        pair_vec_div(S_r.data(),tmp2.data(),R_.data(),comp);
                    }
                    else 
                    {
                        pair_vec_power(S_r.data(),tmp1.data(),j-r,comp);
                        pair_vec_multiplies(tmp2.data(),S_r.data(),S_j.data(),comp);
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_div(S_r.data(),tmp2.data(),tmp1.data(),comp);
                    }
                }
                else
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    O.data()=std::move(S_r.data());
                    return void();
                }
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                if (j-r<2)
                {
                    pair_vec_multiplies(tmp1.data(),R_.data(),tmp2.data(),comp);    
                    prem(tmp2,S_j_1,S_j,v);
                    pair_vec_div(S_r_1.data(),tmp2.data(),tmp1.data(),comp);    
                }        
                else
                {
                    pair_vec_multiplies(R_.data(),tmp2.data(),tmp1.data(),comp);
                    prem(tmp2,S_j_1,S_j,v);
                    pair_vec_div(S_r_1.data(),tmp2.data(),R_.data(),comp);
                }
                if ((j-r) & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=r-1;
            }
            else
            {
                prem(tmp1,S_j_1,S_j,v);
                // std::cout<<"tmp1:"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                // std::cout<<"R_:"<<R_<<std::endl;
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                // std::cout<<"tmp2:"<<tmp2<<std::endl;
                pair_vec_div(S_j_1.data(),tmp1.data(),tmp2.data(),comp);
                --j;
                //  std::cout<<j<<":"<<S_j_1<<std::endl;
                swap(S_j.data(),S_j_1.data());
            }
            
        }
        O.data()=std::move(S_j.data());
        
    }

     
/**
 * @brief  计算f和g的子结式链 要求degree(f,x)>=degree(g,x)>0 ，返回下标与子结式链下标一致
 *  
 * @tparam Tc 
 * @param F 
 * @param G 
 * @param v 
 * @param is_list 
 * @return std::vector<polynomial_<Tc,univariate_priority_order>> 
 */
    template <class Tc>
    std::vector<polynomial_<Tc,univariate_priority_order>> subresultant
        (   
            const polynomial_<Tc,univariate_priority_order>&F,
            const polynomial_<Tc,univariate_priority_order>&G,
            variable v,
            bool is_list=true 
        )
    {
        polynomial_<Tc,univariate_priority_order> O;
        std::vector<polynomial_<Tc,univariate_priority_order>>  L;
        O.clear();
        const univariate_priority_order &comp=G.comp();
       
        int64_t m=get_up_deg(F);
        int64_t l=get_up_deg(G);
        assert(m>=l && l>0);
        
        polynomial_<Tc,univariate_priority_order> S_j_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_j(&comp);
        polynomial_<Tc,univariate_priority_order> S_r_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_r(&comp);
        
        polynomial_<Tc,univariate_priority_order> R_(&comp);
        polynomial_<Tc,univariate_priority_order> tmp1(&comp);
        polynomial_<Tc,univariate_priority_order> tmp2(&comp);
        int64_t j,r,u;
        // std::cout<<"F="<<F<<std::endl;
        // std::cout<<"G="<<G<<std::endl;
        if (l<m)
            u=m;
        else
            u=m+1;
        if (is_list)
        {
            L.resize(u+1,tmp1);
            L[u]=F;
            L[u-1]=G;
        }
        else
            L.resize(1,tmp1);
        if (l<m)
        {
            j=m-1;
            if(j==l)
            {
                prem(S_j,F,G,v);
                --j;
                //  std::cout<<j<<":"<<S_j<<std::endl;
                S_j_1=G;
            }
            else{
                leadcoeff(R_,G);
                if (j-l<2)
                    pair_vec_multiplies(S_r.data(),R_.data(),G.data(),comp);
                else 
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-l,comp);
                    pair_vec_multiplies(S_r.data(),tmp1.data(),G.data(),comp);
                }
                //  std::cout<<l<<":"<<S_r<<std::endl;
                prem(S_r_1,F,G,v);
                if (j-l & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<l-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=l-1;
            }

        }else
        {
            j=m;
            prem(S_j,F,G,v);
            pair_vec_negate(S_j.data());
            //  std::cout<<j-1<<":"<<S_j<<std::endl;
            if (!(--j))
            {
                L[0]=std::move(S_j);
                return L;
            }
            if (S_j.empty())
            {
                if (is_list)
                {
                    L[m-1]=std::move(S_j);
                }
                return L;
            }
            r=get_up_deg(S_j);
            if (r<j)
            {
                
                if (r!=1)
                {
                    leadcoeff(R_,S_j);
                    if (j-r<2)
                        pair_vec_multiplies(S_r.data(),R_.data(),S_j.data(),comp);
                    else 
                    {
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_multiplies(S_r.data(),tmp1.data(),S_j.data(),comp);
                    }
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    if (is_list)
                    {
                        L[m-1]=std::move(S_j);
                    }
                    L[0]=std::move(S_r);
                    return L;
                }
                prem(tmp1,G,S_j,v);
                //std::cout<<"G"<<":"<<G<<std::endl;
                //std::cout<<"S_"<<j<<":"<<S_j<<std::endl;
                
                //  std::cout<<r-1<<":"<<tmp1<<std::endl;
                leadcoeff(R_,G);
                pair_vec_div(S_r_1.data(),tmp1.data(),R_.data(),comp);
                if ((j-r) & 1==1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                j=r-1;
                if (is_list)
                {
                    L[m-1]=std::move(S_j);
                }
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
            }
            else
            {
                prem(tmp1,G,S_j,v);
                leadcoeff(R_,G);
                pair_vec_div(S_j_1.data(),tmp1.data(),R_.data(),comp);
                //  std::cout<<j-1<<":"<<S_j_1<<std::endl;
                --j;
                swap(S_j.data(),S_j_1.data());
            }
        }
        while (j>0)
        {
            if (S_j.empty())
            {
                if (is_list)
                {
                    L[j+1]=std::move(S_j_1);
                    L[j]=std::move(S_j);
                }
                return L;
            }
            r=get_up_deg(S_j);
            if (r<j)
            {
                
                //std::cout<<"R_j"<<":"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                //std::cout<<"R_j+1"<<":"<<R_<<std::endl;
                if (r!=1)
                {
                    leadcoeff(tmp1,S_j);
                    if (j-r<2)
                    {
                        pair_vec_multiplies(tmp2.data(),tmp1.data(),S_j.data(),comp);
                        pair_vec_div(S_r.data(),tmp2.data(),R_.data(),comp);
                    }
                    else 
                    {
                        pair_vec_power(S_r.data(),tmp1.data(),j-r,comp);
                        pair_vec_multiplies(tmp2.data(),S_r.data(),S_j.data(),comp);
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_div(S_r.data(),tmp2.data(),tmp1.data(),comp);
                    }
                }
                else
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    if (is_list)
                    {
                        L[j+1]=std::move(S_j_1);
                        L[j]=std::move(S_j);
                    }
                    L[0]=std::move(S_r);
                    return L;
                }
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                if (j-r<2)
                {
                    pair_vec_multiplies(tmp1.data(),R_.data(),tmp2.data(),comp);    
                    prem(tmp2,S_j_1,S_j,v);
                    pair_vec_div(S_r_1.data(),tmp2.data(),tmp1.data(),comp);    
                }        
                else
                {
                    pair_vec_multiplies(R_.data(),tmp2.data(),tmp1.data(),comp);
                    prem(tmp2,S_j_1,S_j,v);
                    pair_vec_div(S_r_1.data(),tmp2.data(),R_.data(),comp);
                }
                if ((j-r) & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                if (is_list)
                {
                    L[j+1]=std::move(S_j_1);
                    L[j]=std::move(S_j);
                    
                }
                j=r-1;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                
            }
            else
            {
                prem(tmp1,S_j_1,S_j,v);
                // std::cout<<"tmp1:"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                // std::cout<<"R_:"<<R_<<std::endl;
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                // std::cout<<"tmp2:"<<tmp2<<std::endl;
                if (is_list)
                {
                    L[j+1]=std::move(S_j_1);
                }
                pair_vec_div(S_j_1.data(),tmp1.data(),tmp2.data(),comp);
                --j;
                //  std::cout<<j<<":"<<S_j_1<<std::endl;
              
                swap(S_j.data(),S_j_1.data());
            }
            
        }
        if (is_list)
        {
            L[j+1]=std::move(S_j_1);
            L[j]=std::move(S_j);
        }
        else{
            L[0]=std::move(S_j);
        }
        return L;   
    }


    template <class Tc,class comp>
    inline std::vector<polynomial_<Tc,comp>> subresultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v,bool is_list=true)
    {
        assert(comp_consistent(G.comp(),F.comp()));
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        auto O1=subresultant(G1,F1,v,is_list);
        std::vector<polynomial_<Tc,comp>> O(O1.size(),polynomial_<Tc,comp>(G.comp_ptr()));
        for (auto i=0;i<O.size();++i)
        {
            poly_convert(std::move(O1[i]),O[i]);
        }
        return O;
    }

    template <class Tc,class comp>
    polynomial_<Tc,comp> discriminant(const polynomial_<Tc,comp> &G,const variable & v)
    {
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(G,G1);
        auto m=get_up_deg(G1);
        if (m==0)
            return O;
        auto F1=derivative(G1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        resultant(O1,G1,F1,v);
        
        if ((m*(m-1)) & 3)
            O1=-O1;
        O1=O1/leadcoeff(G1);
        poly_convert(std::move(O1),O);
        return O;
    }


    template<class Tc>
    std::vector<polynomial_<Tc,univariate_priority_order>> BezoutMatrix
    (const polynomial_<Tc,univariate_priority_order> &F,
    const polynomial_<Tc,univariate_priority_order> &G)
    {
        auto &xcomp=F.comp();
        // std::cout<<"F="<<F<<";"<<std::endl;
        // std::cout<<"G="<<G<<";"<<std::endl;
        int64_t n=get_up_deg(F);
        int64_t m=get_up_deg(G);
        // std::cout<<"n="<<n<<std::endl;
        assert(n>=m);
        polynomial_<Tc,univariate_priority_order> tmp(&xcomp),f1(&xcomp),f2(&xcomp),f3(&xcomp),f4(&xcomp);
        auto coF=coeff(F);
        auto coG=coeff(G);
        
    
        // std::cout<<"F:"<<coF<<std::endl;
        // std::cout<<"G:"<<coG<<std::endl;
        
        std::vector<polynomial_<Tc,univariate_priority_order>>  
            B(n*n,tmp);
        for(auto i=1;i<m+1;i++)
        {
            f1.clear();f2.clear();f3.clear();f4.clear();
            for(auto r=i;r<m+1;r++) 
            {
                tmp.clear();
                for(auto &t: coG[m-r])
                    tmp.push_back({__change_up_monomial_var_deg(t.first,n-r),t.second});
                f1=f1+tmp;
            }
            // f1=f1+coG[m-r]*pow(x,n-r);
            for(auto r=i;r<n+1;r++) 
            {
                tmp.clear();
                for(auto &t: coF[n-r])
                    tmp.push_back({__change_up_monomial_var_deg(t.first,n-r),t.second});
                f3=f3+tmp;
            }
            // f3=f3+coF[n-r]*pow(x,n-r);
            for(auto r=1;r<i+1;r++)
            {
                tmp.clear();
                for(auto &t: coF[n-r+1])
                    tmp.push_back({__change_up_monomial_var_deg(t.first,i-r),t.second});
                f2=f2+tmp;
                // f2=f2+coF[n-r+1]*pow(x,i-r);
                tmp.clear();
                for(auto &t: coG[m-r+1])
                    tmp.push_back({__change_up_monomial_var_deg(t.first,i-r),t.second});
                f4=f4+tmp;
                // f4=f4+coG[m-r+1]*pow(x,i-r);
            }
            tmp=f1*f2-f3*f4;
            
            auto cot=coeff(tmp);
            for(auto j=1;j<n+1;j++) if (n-j<cot.size()) B[(m-i)*n+j-1]=cot[n-j];
        }
        //case 2
        for(auto i=m+1;i<n+1;i++)
            for(auto j=1;j<m+1+1;j++)
                B[(i-1)*n+(i-m-1+j-1)]=coG[m-j+1];
        return B;

    }


    template <class Tc>
    polynomial_<Tc,univariate_priority_order> det(std::vector<polynomial_<Tc,univariate_priority_order>> & A,size_t n)
    {
        //求行列式值
        
        assert(A.size()>=n*n && n>0);
        if (n == 1)
            return A[0];
        std::vector<polynomial_<Tc,univariate_priority_order>>B((n-1)*(n-1));//创建n-1阶的代数余子式阵B  
        int mov = 0;//判断行是否移动   
        polynomial_<Tc,univariate_priority_order> sum(A[0].comp_ptr()) ;//sum为行列式的值  
        for (int arow = 0; arow<n; arow++) // A的行数把矩阵A(nn)赋值到A(n-1,n-1)  
        {
            if (A[arow*n].empty()) continue;
            for (int brow = 0; brow<n - 1; brow++)//把A阵第一列各元素的代数余子式存到bb  
            {    
                mov = arow > brow ? 0 : 1; //B中小于arow的行，同行赋值，等于的错过，大于的加一  
                for (int j = 0; j<n - 1; j++)  //从A的第二列赋值到第n列  
                {
                    B[brow*(n-1)+j] = A[(brow + mov)*n+ j + 1];
                }
            }
            sum = (arow % 2 == 0 ? sum+A[arow*n] * det(B,n-1): sum-A[arow*n] * det(B,n-1)); //因为列数为0，所以行数是偶数时候，代数余子式为1.  
        }
        return sum;
    }



    
}
#endif