/**
 * @file upolynomial.cc
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义upolynomial
*/
#include <clpoly/upolynomial.hh>
namespace clpoly{
    uless uless::init;
    void poly_convert(const upolynomial_<ZZ>& p_in,upolynomial_<QQ> & p_out)
    {
        p_out.clear();
        p_out.reserve(p_in.size());
        for(auto &i:p_in)
        {
            p_out.push_back({i.first,i.second});
        }
        p_out.normalization();
    }
    void poly_convert(const upolynomial_<QQ>& p_in,upolynomial_<ZZ> & p_out)
    {
        p_out.clear();
        p_out.reserve(p_in.size());
        if (p_in.empty())
            return void();

        ZZ den=1;
        for (auto &i:p_in)
        {
            den=lcm(den,i.second.get_den());
        }
        for(auto &i:p_in)
        {
            p_out.push_back({i.first,i.second.get_num()*(den/i.second.get_den())});
        }
        p_out.normalization();
    }
}