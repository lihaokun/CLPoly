/**
 * @file random.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 
 * 
 */

#ifndef CLPOLY_RANDOM_HH
#define CLPOLY_RANDOM_HH

#include <clpoly/polynomial.hh>
#include <clpoly/associatedgraph.hh>
#include <list>
#include <string>
#include <random>
#include <algorithm>


namespace clpoly{

    // ---- combinatorial helpers for monomial unranking ----
    namespace detail{
        // C(n, k)
        inline uint64_t binomial(uint64_t n, uint64_t k) {
            if (k > n) return 0;
            if (k == 0 || k == n) return 1;
            if (k > n - k) k = n - k;
            uint64_t r = 1;
            for (uint64_t i = 0; i < k; ++i)
                r = r * (n - i) / (i + 1);
            return r;
        }

        // Number of monomials in n_vars variables with total degree <= max_deg
        // = C(max_deg + n_vars, n_vars)
        inline uint64_t monomial_count_upto(uint64_t n_vars, uint64_t max_deg) {
            return binomial(max_deg + n_vars, n_vars);
        }

        // Number of monomials in n_vars variables with total degree == deg
        // = C(deg + n_vars - 1, n_vars - 1)
        inline uint64_t monomial_count_exact(uint64_t n_vars, uint64_t deg) {
            return binomial(deg + n_vars - 1, n_vars - 1);
        }

        // Unrank index in [0, monomial_count_upto(n_vars, max_deg)-1]
        // to an exponent vector of a monomial with degree <= max_deg.
        // Ordering: ascending first-variable exponent, then recursive.
        inline std::vector<int64_t> unrank_monomial_upto(uint64_t index, uint64_t n_vars, uint64_t max_deg) {
            std::vector<int64_t> exp(n_vars, 0);
            uint64_t rem = max_deg;
            for (uint64_t i = 0; i < n_vars; ++i) {
                for (uint64_t e = 0; e <= rem; ++e) {
                    uint64_t cnt = (i + 1 < n_vars)
                        ? monomial_count_upto(n_vars - i - 1, rem - e)
                        : 1;
                    if (index < cnt) {
                        exp[i] = static_cast<int64_t>(e);
                        rem -= e;
                        break;
                    }
                    index -= cnt;
                }
            }
            return exp;
        }

        // Unrank index in [0, monomial_count_exact(n_vars, deg)-1]
        // to an exponent vector of a monomial with degree == deg.
        inline std::vector<int64_t> unrank_monomial_exact(uint64_t index, uint64_t n_vars, uint64_t deg) {
            std::vector<int64_t> exp(n_vars, 0);
            uint64_t rem = deg;
            for (uint64_t i = 0; i + 1 < n_vars; ++i) {
                for (uint64_t e = 0; e <= rem; ++e) {
                    uint64_t cnt = monomial_count_exact(n_vars - i - 1, rem - e);
                    if (index < cnt) {
                        exp[i] = static_cast<int64_t>(e);
                        rem -= e;
                        break;
                    }
                    index -= cnt;
                }
            }
            exp[n_vars - 1] = static_cast<int64_t>(rem);
            return exp;
        }
    } // namespace detail
    inline std::vector<size_t> random_select(size_t m,size_t n,bool no_repeat=true) // 0,..,m-1 中选 n 个
    {
        std::vector<size_t> v;
        if (no_repeat && n>=m)
        {
            n=m;
            for(size_t i=0;i!=n;++i)
                v.push_back(i);
            return v;
        
        }
        std::random_device rd; 
        std::mt19937 gen(rd());
        size_t tmp=no_repeat?m-n:m-1;
        std::uniform_int_distribution<> dis(0, tmp);
        for (size_t i=0;i!=n;++i)
        {
            v.push_back(dis(gen));
        }
        std::sort(v.begin(),v.end());
        if (no_repeat)
        {
            for (size_t i=0;i!=n;++i)
            {
                v[i]+=i;
            }
        }
        return v;
    }
    
    /**
     * Generate a random sparse polynomial using combinatorial unranking.
     *
     * Algorithm (based on REDUCE's randpoly):
     *   1. Compute the monomial space: C(deg+n, n) monomials with degree <= deg.
     *   2. Pick one monomial of degree exactly `deg` (guarantees polynomial degree).
     *   3. Pick `len-1` distinct monomials from the remaining space via unranking.
     *   4. Assign random coefficients from [coeff.first, coeff.second].
     *
     * All monomials are guaranteed unique — no rejection sampling needed.
     */
    template<class Tc>
    polynomial_<Tc> random_polynomial(const std::vector<variable> & v,uint64_t deg,uint64_t len,std::pair<int,int>  coeff,bool is_add_num=false)
    {
        if (coeff.first > coeff.second) std::swap(coeff.first,coeff.second);
        uint64_t n = v.size();
        uint64_t total = detail::monomial_count_upto(n, deg);
        uint64_t hi    = detail::monomial_count_exact(n, deg); // degree == deg
        uint64_t lo    = total - hi;                           // degree <  deg

        if (len > total) len = total;
        if (len == 0) return polynomial_<Tc>();

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(coeff.first, coeff.second);

        // Generate a non-zero coefficient
        auto nonzero_coeff = [&]() -> Tc {
            Tc c(dis(gen));
            while (!c) c = Tc(dis(gen));
            return c;
        };

        // Helper: exponent vector -> monomial
        std::vector<std::pair<variable,int64_t>> m_tmpl;
        for (auto &i:v) m_tmpl.emplace_back(i, 0);
        auto make_mono = [&](const std::vector<int64_t>& exp) {
            auto m = m_tmpl;
            for (size_t i = 0; i < n; ++i) m[i].second = exp[i];
            return monomial(m);
        };

        std::vector<std::pair<monomial,Tc>> p1;
        p1.reserve(len + (is_add_num ? 1 : 0));

        // Step 1: one monomial of degree exactly `deg`
        std::uniform_int_distribution<uint64_t> hi_dis(0, hi - 1);
        uint64_t first_idx = hi_dis(gen);
        p1.emplace_back(
            make_mono(detail::unrank_monomial_exact(first_idx, n, deg)),
            nonzero_coeff());

        // Step 2: len-1 distinct monomials from the remaining total-1 slots
        if (len > 1) {
            auto indices = random_select(total - 1, len - 1, true);
            for (auto idx : indices) {
                std::vector<int64_t> exp;
                if (idx < lo) {
                    exp = detail::unrank_monomial_upto(idx, n, deg - 1);
                } else {
                    uint64_t exact_idx = idx - lo;
                    if (exact_idx >= first_idx) ++exact_idx;
                    exp = detail::unrank_monomial_exact(exact_idx, n, deg);
                }
                p1.emplace_back(make_mono(exp), nonzero_coeff());
            }
        }

        // Step 3: add constant term if requested (only if not already present)
        if (is_add_num) {
            bool has_const = false;
            for (auto &t : p1)
                if (t.first.empty()) { has_const = true; break; }
            if (!has_const)
                p1.push_back({{}, nonzero_coeff()});
        }
        return polynomial_<Tc>(p1);
    }

    template<class Tc>
    polynomial_<Tc> random_polynomial(const std::vector<variable> & v,uint64_t deg,double p,int up,int down)
    {
        if (down > up) std::swap(up,down);
        std::vector<std::pair<monomial,Tc>> p1;
        std::random_device rd; 
        std::mt19937 gen(rd());
        auto size=v.size();
        std::bernoulli_distribution mp(p);
        std::uniform_int_distribution<> dis(down, up);
        std::vector<std::pair<variable,int64_t>>  m;
        for (auto &i:v)
            m.emplace_back(i,0);
        for(int d=deg;d>=0;--d)
        {
            
            if (d)
            {
                int64_t i0=0;
                int64_t sum=d;
                m.begin()->second=d;
                while(1)
                {
                    if (mp(gen))
                    {
                        p1.emplace_back(monomial(m),Tc(dis(gen)));
                    }
                    while (i0>=0)
                    {
                        if (m[i0].second && i0<size-1)
                        {
                            --m[i0].second;
                            ++i0;
                            m[i0].second=d-sum+1;
                            sum=d;
                            break;
                        }
                        else
                        {
                            sum-=m[i0].second;
                            m[i0].second=0;
                            --i0;
                        }
                    }   
                    if (i0<0)   
                        break;

                }
            }
            else
            {
                if (mp(gen))
                {
                    p1.push_back({{},Tc(dis(gen))});
                }
            }
            
        }
        return polynomial_<Tc>(p1);
    }

    template <class T>
    std::vector<T> RandomSample(const  std::vector<T> & l,uint64_t n)
    {
        
        std::random_device rd; 
        std::mt19937 gen(rd());
        auto size=l.size();
        std::vector<T>  lout;
        lout.reserve(std::min(n,size));
        std::vector<bool> bl;
        bl.reserve(size);
        uint64_t select;
        for (uint64_t i=size;i>0;--i)
            bl.push_back(true);
        for (uint64_t i=std::min(n,size);i>0;--i)
        {
            if (i==1)
            {
                select = 1;
            }
            else{
                std::uniform_int_distribution<uint64_t> dis(1, i);
                select=dis(gen);
            }
            for (uint64_t j=0;j<size;++j)
            {
                if (bl[j])
                {
                    --select;
                    if (!select)
                    {
                        bl[j]=false;
                        lout.push_back(l[j]);
                        break;
                    }
                }
            }
        }
        return lout;
    }

    template <class node>
    graph<node> random_graph(const std::vector<node> & nodes,double p)
    {
        graph<node> G;
        for (auto & i:nodes)
            G.add_node(i);
        size_t n=G.size();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::bernoulli_distribution d(p);
        for (size_t i=0;i<n;++i)
        {
            for (size_t j=i+1;j<n;++j)
                if (d(gen))
                    G.add_edge_index(i,j);
        }
        return G;
    }
    /**
     * random_polynomials<Tc>
     * 生成一组多项式
     * 变量相关图的稀疏性是 p1
     * 项稀疏性是 p2^std::max(n-1,1) n为这个多项式变量数
    **/
    template<class Tc>
    std::vector<polynomial_<Tc>> random_polynomials(const std::vector<variable> & vars,uint64_t deg,double p1,double p2,int up,int down)
    {
        auto G=random_graph<variable>(vars,p1);
        // std::cout<<G<<std::endl;
        std::vector<polynomial_<Tc>> l;
        std::vector<variable> v_l;
        
        for (auto &v:G.nodes())
        {
            v_l.clear();
            v_l.push_back(v);
            auto i=G.index(v);
            double p=1;
            for (auto &j:G.adjacency_list()[i])
            {
                if (j!=i)
                {
                    bool b=true;
                    auto node_j=G.node(j);
                    for (auto &k:v_l)
                    {
                        if (!G.is_edge(k,node_j))
                        {
                            b=false;
                            break;
                        }
                    }
                    if (b)
                    {
                        v_l.push_back(node_j);
                        p*=p2;
                    }
                }
            }
            if(v_l.size()>1)
                l.push_back(random_polynomial<clpoly::ZZ>(v_l,deg,p,up,down));
        }
        return l;
    }
}
#endif