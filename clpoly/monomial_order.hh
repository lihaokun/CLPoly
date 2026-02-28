/**
 * @file monomial_order.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义类：单项序和单项压缩算法
 * 
 */
#ifndef CLPOLY_MONOMIAL_ORDER_HH
#define CLPOLY_MONOMIAL_ORDER_HH
#include <functional>
#include <cassert>
#include <map>
#include <clpoly/basic_monomial.hh>
namespace clpoly{
    
    class less
    {
        public:
            constexpr bool operator()(const variable & v1,const variable & v2) const 
            {
                return  v1<v2; 
            }
            constexpr bool operator==(const less &g1)const {return true;}
            constexpr bool operator!=(const less &g1)const {return false;}  
    };
    class custom_var_order
    {
        public:
            std::map<variable,int64_t> v_map;
            custom_var_order(){
                // v_map[variable()]=0;
                }
            custom_var_order(const std::vector<variable> &v)
            {
                // v_map[variable()]=0;
                for (int64_t i=0;i<v.size();++i)
                    v_map[v[i]]=i+1;
            }
            // 允许用初始化列表构造
            custom_var_order(std::initializer_list<variable> lst)
                : custom_var_order(std::vector<variable>(lst)) {}
            int64_t order(variable v) const
            {
                // if (v.serial()==0) return 0;
                auto ptr=v_map.find(v);
                if (ptr==v_map.end())
                    return v_map.size()+1;
                else
                    return ptr->second;
            }
            inline bool operator()(const variable & v1,const variable & v2) const 
            {
                auto v1_= this->order(v1);
                auto v2_=this->order(v2);
                return  (v1_<v2_ || (v1_==v2_ && v1<v2)); 
            }
            inline bool operator==(const custom_var_order &g1)const {return v_map==g1.v_map;}
            inline bool operator!=(const custom_var_order &g1)const {return v_map!=g1.v_map;}
    };
    class univariate_order
    {
        public:     
            variable v;
            univariate_order(){}
            univariate_order(variable _v):v(_v){}
            constexpr bool operator()(const variable & v1,const variable & v2) const 
            {
                return  (v2!=v && (v1==v || v1<v2)); 
            }
            constexpr bool operator==(const univariate_order &g1)const {return v==g1.v;}
            constexpr bool operator!=(const univariate_order &g1)const {return v!=g1.v;}
    };

    template<class var_comp>
    class lex_
    {
        public:
            var_comp comp;
            static lex_ init;
            lex_():comp(){}
            lex_(var_comp c):comp(std::move(c)){}
            constexpr bool operator()(const variable & v1,const variable & v2) const {return comp(v1,v2);}
            inline bool operator()(const basic_monomial<lex_> &m1,const basic_monomial<lex_> &m2)const {return pair_vec_comp(m1.data(),m2.data(),*this);}
            inline bool operator==(const lex_ &g1)const {return comp==g1.comp;}
            inline bool operator!=(const lex_ &g1)const {return comp!=g1.comp;}
    };
     template<class var_comp>
    lex_<var_comp> lex_<var_comp>::init;
    typedef lex_<less> lex;
    // struct lex
    // {
    //     constexpr bool operator()(const variable & v1,const variable & v2) const {return v1<v2;}
    //     //constexpr bool operator()(uint64_t v1,uint64_t v2) const {return v1<v2;}
    //     inline bool operator()(const basic_monomial<lex> &m1,const basic_monomial<lex> &m2)const {return pair_vec_comp(m1.data(),m2.data(),m1.comp());}
    //     constexpr bool operator==(const lex &g1)const {return true;}
    //     constexpr bool operator!=(const lex &g1)const {return false;}
        
    // };
    template<class var_comp>
    class grlex_
    {
        public:
            var_comp comp;
            static grlex_ init;
            grlex_():comp(){}
            grlex_(var_comp c):comp(std::move(c)){}
            // template<class T>
            // grlex_(const T& c):comp(c){}
            
            constexpr bool operator()(const variable & v1,const variable & v2) const {return comp(v1,v2);}
            inline bool operator()(const basic_monomial<grlex_> &m1,const basic_monomial<grlex_> &m2)const {return m1.deg()>m2.deg() || (m1.deg()==m2.deg() && pair_vec_comp(m1.data(),m2.data(),*this));}
            inline bool operator==(const grlex_ &g1)const {return comp==g1.comp;}
            inline bool operator!=(const grlex_ &g1)const {return comp!=g1.comp;}
    };
    template<class var_comp>
    grlex_<var_comp> grlex_<var_comp>::init;
    typedef grlex_<less> grlex;
    // struct grlex
    // {
    //     constexpr bool operator()(const variable & v1,const variable & v2) const {return v1<v2;}
    //     //constexpr bool operator()(uint64_t v1,uint64_t v2) const {return v1>v2;}
    //     inline bool operator()(const basic_monomial<grlex> &m1,const basic_monomial<grlex> &m2)const {return m1>m2;}
    //     constexpr bool operator==(const grlex &g1)const {return true;}
    //     constexpr bool operator!=(const grlex &g1)const {return false;}
    // };
    // template<class var_comp>
    // struct grevlex_
    // {
    //     var_comp comp;
    //     grevlex_():comp(){}
    //     grevlex_(var_comp c):comp(std::move(c)){}
    //     constexpr bool operator()(const variable & v1,const variable & v2) const {return comp(v1,v2);}
    //     inline bool operator()(const basic_monomial<grlex_> &m1,const basic_monomial<grlex_> &m2)const {return m1.deg()>m2.deg() || (m1.deg()==m2.deg() && pair_vec_comp(m1.data(),m2.data(),*this));}
    //     inline bool operator==(const grlex_ &g1)const {return comp==g1.comp;}
    //     inline bool operator!=(const grlex_ &g1)const {return comp!=g1.comp;}
    // };
    // struct grevlex
    // {
    //     constexpr bool operator()(const variable & v1,const variable & v2) const {return v1>v2;}
    //     inline bool operator()(const basic_monomial<grevlex> &m1,const basic_monomial<grevlex> &m2)const {return m1>m2;}
    //     constexpr bool operator==(const grevlex &g1)const {return true;}
    //     constexpr bool operator!=(const grevlex &g1)const {return false;}
    // };
    // template <class comp>
    // constexpr const int64_t get_up_deg(const basic_monomial<comp>& m1){return m1.empty()?0:m1.back().second;}
    

    
    
    // struct univariate_priority_order;
    // constexpr const int64_t get_up_deg(const basic_monomial<univariate_priority_order>& m1);

    // struct univariate_priority_order
    // {
    //     variable v;
    //     static univariate_priority_order init;
    //     univariate_priority_order():v(){}
    //     univariate_priority_order(variable _v):v(_v){}
    //     constexpr bool operator()(const variable & v1,const variable & v2) const {return(v1!=v && (v2==v || v1<v2));}
    //     inline bool operator()(const basic_monomial<univariate_priority_order> &m1,const basic_monomial<univariate_priority_order> &m2)const 
    //     {
    //         assert(m1.comp().v==m2.comp().v);
    //         auto d1=get_up_deg(m1);
    //         auto d2=get_up_deg(m2);
    //         return (d1>d2 || (d1==d2 && pair_vec_comp(m1.data(),m2.data(),m1.comp())));
    //     }
    //     constexpr bool operator==(const univariate_priority_order &g1)const {return v==g1.v;}
    //     constexpr bool operator!=(const univariate_priority_order &g1)const {return v!=g1.v;}
    //     constexpr const variable & var() const{return v;}
    // };
    // univariate_priority_order univariate_priority_order::init;

    typedef lex_<univariate_order> univariate_priority_order;
    inline  variable get_up_var(const univariate_priority_order & comp)
    {return comp.comp.v;}
    constexpr int64_t get_up_deg(const basic_monomial<univariate_priority_order>& m1)
    {return (m1.empty() || m1.front().first!=get_up_var(m1.comp()))?0:m1.front().second;}


/*compression*/

    template<class Tc1,class Tc2,class compare,class compare2>
    constexpr bool __is_monomial_can_compression(
        const std::vector<std::pair<basic_monomial<compare>,Tc1>> & v1_,
        const std::vector<std::pair<basic_monomial<compare>,Tc2>> & v2_,
        const compare2 & comp,
        std::list<variable>& vars,
        int delta=0
    )
    {
        return false;
    }
    template<class compare>
    constexpr uint64_t __monomial_compression(const basic_monomial<compare> & m,const std::list<variable>& vars)
    {
        return 0;
    }
    template<class compare>
    constexpr void __monomial_decompression(uint64_t mc,basic_monomial<compare> & m,const std::list<variable>& vars,const compare * comp)
    {}
    template<class compare>
    constexpr uint64_t __monomial_compression_div_mold(const compare * m,size_t vars_size)
    {
        return 0;
    }
     
    /*grlex*/
    template<class Tc1,class Tc2, class var_comp>
    bool __is_monomial_can_compression(
        const std::vector<std::pair<basic_monomial<grlex_<var_comp>>,Tc1>> & v1_,
        const std::vector<std::pair<basic_monomial<grlex_<var_comp>>,Tc2>> & v2_,
        const grlex_<var_comp> & comp,
        std::list<variable>& vars,
        int delta=0
    )
    {
        for (auto & i:v1_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
        for (auto & i:v2_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
        vars.clear();
        __pair_vec_variables(v1_,vars);
        __pair_vec_variables(v2_,vars);
        if (delta)
            if (vars.size()>1 && std::max(v1_.begin()->first.deg(),v2_.begin()->first.deg())>=(uint64_t(1)<<(64/(vars.size()+1)-1))) return false;
        else
            if (vars.size()>1 && v1_.begin()->first.deg()+v2_.begin()->first.deg()>=(uint64_t(1)<<(64/(vars.size()+1)))) return false;
        return true;
    }
    template<class var_comp>
    inline uint64_t __monomial_compression(const basic_monomial<grlex_<var_comp>> & m,const std::list<variable>& vars)
    {
        uint64_t mc=0;
        if (m.empty() || vars.empty() )return mc;
        if (vars.size()>1)
        {
            uint l=64/(vars.size()+1);
            mc=m.deg();
            auto m_ptr=m.begin();
            for (auto &i:vars)
            {
                //if (m_ptr==m.end())

                mc<<=l;
                if (m_ptr!=m.end()  && i == m_ptr->first)
                {
                    mc+=m_ptr->second;
                    ++m_ptr;
                }
            }
        }
        else
        {
            mc=m.deg();
        }
        return mc;
    }
    template<class var_comp>
    inline void __monomial_decompression(uint64_t mc,basic_monomial<grlex_<var_comp>> & m,const std::list<variable>& vars,const grlex_<var_comp> * comp_ptr)
    {
        m.clear();
        m.comp(comp_ptr);
        if (!mc || vars.empty()) return void();
        m.reserve(vars.size());
        if (vars.size()>1)
        {
            uint l=64/(vars.size()+1);
            uint ll=l*vars.size();
            uint64_t mod=(uint64_t(1)<<(l*(vars.size()+1)))-(uint64_t(1)<<ll);
            m.deg()=(mc & mod)>>ll;
            uint64_t part;
            for (auto &i:vars)
            {
                mc<<=l;
                part=(mc & mod)>>ll;
                if (part!=0)
                    m.data().push_back({i,part});
            }
        }
        else
        {
            m.push_back({*vars.begin(),mc});
        }
    }
    template<class var_comp>
    inline uint64_t __monomial_compression_div_mold(const grlex_<var_comp>* m,size_t vars_size)
    {
        uint l=64/(vars_size+1);
        uint64_t mold=0;
        for (uint i=0;i<=vars_size;++i)
        {
            mold<<=1;
            mold+=1;
            mold<<=(l-1);
        }
        return mold;
    }




    // /*univariate_priority_order*/
    // template<class Tc1,class Tc2>
    // bool __is_monomial_can_compression(
    //     const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc1>> & v1_,
    //     const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc2>> & v2_,
    //     const univariate_priority_order & comp,
    //     std::list<variable>& vars,
    //     int delta=0
    // )
    // {
    //     uint64_t deg1=0,deg2=0;
    //     for (auto & i:v1_)
    //         for (auto & j:i.first)
    //             if (j.second<0)
    //                 return false;
    //             else if(j.second>deg1)
    //                 deg1=j.second;

                
    //     for (auto & i:v2_)
    //         for (auto & j:i.first)
    //             if (j.second<0)
    //                 return false;
    //             else if(j.second>deg2)
    //                 deg2=j.second;
    //     vars.clear();
    //     __pair_vec_variables(v1_,vars);
    //     __pair_vec_variables(v2_,vars);
    //     if (delta)
    //         if (vars.size()>1  && std::max(deg1,deg2)>=(uint64_t(1)<<(64/vars.size()-1))) return false;
    //     else    
    //         if (vars.size()>1  && deg1+deg2>=(uint64_t(1)<<(64/vars.size()))) return false;
    //     return true;
    // }
    // template<>
    // inline uint64_t __monomial_compression(const basic_monomial<univariate_priority_order> & m,const std::list<variable>& vars)
    // {
    //     uint64_t mc=0;
    //     if (m.empty() || vars.empty() )return mc;
    //     uint l=64/(vars.size());
    //     variable v=m.comp().v;
    //     auto m_ptr=m.begin();
    //     if (vars.back()==v)
    //     {
    //         mc=get_up_deg(m);
    //     }
    //     for (auto &i:vars)
    //     {
    //         if (i!=v)
    //         {
    //             mc<<=l;
    //             if (m_ptr!=m.end()  && i == m_ptr->first)
    //             {
    //                 mc+=m_ptr->second;
    //                 ++m_ptr;
    //             }
    //         }
    //     }
    //     return mc;
    // }
    // template<>
    // inline void __monomial_decompression(uint64_t mc,basic_monomial<univariate_priority_order> & m,const std::list<variable>& vars,const univariate_priority_order * comp_ptr)
    // {
    //     m.clear();
    //     m.comp(comp_ptr);
    //     if (!mc || vars.empty()) return void();
    //     m.reserve(vars.size());
    //     variable v=comp_ptr->v;
    //     uint l=64/(vars.size());
    //     uint ll=l*(vars.size()-1);
    //     uint64_t mod;
    //     if (l==64)
    //         mod=-1;
    //     else
    //         mod=(uint64_t(1)<<(l*(vars.size())))-(uint64_t(1)<<ll);
    //     uint64_t deg;
    //     if (vars.back()==v)
    //     {
    //         deg=(mc & mod)>>ll;
    //         mc<<=l;
    //     }
    //     uint64_t part;
    //     for (auto &i:vars)
    //     {
    //         if (i!=v)
    //         {
    //             part=(mc & mod)>>ll;
    //             if (part!=0)
    //                 m.push_back({i,part});
    //             mc<<=l;
    //         }
    //         else
    //         {
    //             if (deg !=0)
    //                 m.push_back({v,deg});
    //         }   
    //     }

    // }
    // template<>
    // inline uint64_t __monomial_compression_div_mold(const univariate_priority_order* comp,size_t vars_size)
    // {
    //     uint l=64/(vars_size);
    //     uint64_t mold=0;
    //     for (uint i=0;i<vars_size;++i)
    //     {
    //         mold<<=1;
    //         mold+=1;
    //         mold<<=(l-1);
    //     }
    //     return mold;
    // }




    /*lex*/
    template<class Tc1,class Tc2,class var_order>
    bool __is_monomial_can_compression(
        const std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc1>> & v1_,
        const std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc2>> & v2_,
        const lex_<var_order> & comp,
        std::list<variable>& vars,
        int delta=0
    )
    {
        uint64_t deg1=0,deg2=0;
        for (auto & i:v1_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
                else if(j.second>deg1)
                    deg1=j.second;

                
        for (auto & i:v2_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
                else if(j.second>deg2)
                    deg2=j.second;
        vars.clear();
        __pair_vec_variables(v1_,vars);
        __pair_vec_variables(v2_,vars);
        if (delta)
            if (vars.size()>1  && std::max(deg1,deg2)>=(uint64_t(1)<<(64/vars.size()-1))) return false;
        else    
            if (vars.size()>1  && deg1+deg2>=(uint64_t(1)<<(64/vars.size()))) return false;
        return true;
    }
    template<class var_order>
    inline uint64_t __monomial_compression(const basic_monomial<lex_<var_order>> & m,const std::list<variable>& vars)
    {
        uint64_t mc=0;
        if (m.empty() || vars.empty() )return mc;
        uint l=64/(vars.size());
        auto m_ptr=m.begin();
        for (auto &i:vars)
        {
            mc<<=l;
            if (m_ptr!=m.end()  && i == m_ptr->first)
            {
                mc+=m_ptr->second;
                ++m_ptr;
            }
        }
        return mc;
    }
    template<class var_order>
    inline void __monomial_decompression(uint64_t mc,basic_monomial<lex_<var_order>> & m,const std::list<variable>& vars,const lex_<var_order> * comp_ptr)
    {
        m.clear();
        m.comp(comp_ptr);
        if (!mc || vars.empty()) return void();
        m.reserve(vars.size());
        uint l=64/(vars.size());
        uint ll=l*(vars.size()-1);
        uint64_t mod;
        if (l==64)
            mod=-1;
        else
            mod=(uint64_t(1)<<(l*(vars.size())))-(uint64_t(1)<<ll);
        uint64_t deg;
        uint64_t part;
        for (auto &i:vars)
        {
            part=(mc & mod)>>ll;
            if (part!=0)
                m.push_back({i,part});
            mc<<=l;
          
        }

    }
    template<class var_order>
    inline uint64_t __monomial_compression_div_mold(const lex_<var_order>* comp,size_t vars_size)
    {
        uint l=64/(vars_size);
        uint64_t mold=0;
        for (uint i=0;i<vars_size;++i)
        {
            mold<<=1;
            mold+=1;
            mold<<=(l-1);
        }
        return mold;
    }

}
#endif
