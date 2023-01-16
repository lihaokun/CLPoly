/**
 * @file interval.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义interval class
*/

#include <mpria.h>
#include <sstream>
#include <iostream>
#include <clpoly/number.hh>
#include <clpoly/polynomial_.hh>

#include <clpoly/realroot.hh>
#include <clpoly/interval.hh>
namespace clpoly{
    class g_z{
        public:
        bool operator()(const QQ& i) const
        {
            return i>0;
        }
    };
    
    class ge_z{
        public:
        bool operator()(const QQ& i)const
        {
            return i>=0;
        }
    };
    class l_z{
        public:
        bool operator()(const QQ& i)const
        {
            return i<0;
        }
    };
    class le_z{
        public:
        bool operator()(const QQ& i)const
        {
            return i<=0;
        }
    };
    
    interval feasible_range(const interval_upoly &p,char op)
    {

        std::cout << p<<std::endl;
        upolynomial_<QQ> p_up;
        upolynomial_<QQ> p_down;
        upolynomial_<QQ> n_up;
        upolynomial_<QQ> n_down;
        if (op=='>')
        {
            for (auto &i:p)
            {
                p_up.push_back({i.first,i.second.get_r()});
                if (i.first.deg() % 1 )
                {
                    n_up.push_back({i.first,i.second.get_r()});
                }
                else{
                    n_up.push_back({i.first,-i.second.get_l()});
                }


            }    
            return good_range(p_up,g_z()) || (-good_range(n_up,g_z())); 
        }  
        if (op=='<')
        {
            for (auto &i:p)
            {
                p_down.push_back({i.first,i.second.get_l()});
                if (i.first.deg() % 1 )
                {
                    n_down.push_back({i.first,i.second.get_l()});
                }
                else{
                    n_down.push_back({i.first,-i.second.get_r()});
                }


            }    
            return good_range(p_down,l_z()) || (-good_range(n_down,l_z()));
        }
        if (op=='=')
        {
            for (auto &i:p)
            {
                p_down.push_back({i.first,i.second.get_l()});
                p_up.push_back({i.first,i.second.get_r()});

                if (i.first.deg() % 1 )
                {
                    n_up.push_back({i.first,i.second.get_r()});
                    n_down.push_back({i.first,i.second.get_l()});
                }
                else{
                    n_up.push_back({i.first,-i.second.get_l()});
                    n_down.push_back({i.first,-i.second.get_r()});
                }
            } 
            
            return (good_range(p_up,ge_z())&& good_range(p_down,le_z())) || -(good_range(n_up,ge_z())&& good_range(n_down,le_z())); 
        }
        
        throw std::invalid_argument("feasible_range 不支持的op "+std::string(1,op));
        
    }
    
}
