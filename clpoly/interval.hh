/**
 * @file interval.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义interval class
*/

#ifndef CLPOLY_INTERVAL_HH
#define CLPOLY_INTERVAL_HH
#include <mpria.h>
#include <sstream>
#include <iostream>
#include <clpoly/number.hh>
#include <clpoly/polynomial_.hh>
#include <clpoly/realroot.hh>
namespace clpoly{
    class interval
    {
        private:
            mpri_t _I;

            inline static std::ostream& _to_stream(std::ostream& stream, const mpq_t q)
            {
                if (mpria_mpq_is_nan(q))
                {
                    stream<<"NAN";
                }
                else if (mpria_mpq_is_infinite(q)){
                    if (mpria_mpq_sgn(q)<0)
                        stream<<"-";
                    stream<<"INF";
                }
                else{
                    stream<<q;
                }
                return stream;
            }
        public:

            inline interval(){
                mpri_init(_I);
                // MPRIA_MPQ_SET_NEG_INF(mpri_lepref(_I));
                // MPRIA_MPQ_SET_POS_INF(mpri_repref(_I));
            }

            inline interval(const interval& I){
                mpri_init(_I);
                mpri_set(this->_I,I._I);
            }
            inline interval(interval&& I){
                mpri_init(_I);
                mpri_swap(this->_I,I._I);
            }
            inline interval(const int& i)
            {
                mpri_init(_I);
                mpq_set_si(mpri_lepref(this->_I),i,1);
                mpq_set_si(mpri_repref(this->_I),i,1);
            }
            inline interval(const ZZ& z)
            {
                mpri_init(_I);
                mpq_set_z(mpri_lepref(this->_I),z.get_mpz_t());
                mpq_set_z(mpri_repref(this->_I),z.get_mpz_t());
            }
            
            inline interval(const QQ& q){
                mpri_init(_I);
                MPRI_SET_Q(_I,q.get_mpq_t());
            }
            inline interval(const QQ& l,const QQ& r){
                mpri_init(_I);
                if (QQ_cmp(l,r)>0)
                {
                    MPRI_SET_NAN(_I)
                    return;
                }
                this->set_l(l);
                this->set_r(r);
                // MPRI_SET_Q(_I,q.get_mpq_t());
            }

            inline  virtual ~interval(){
                mpri_clear(_I);
            }
            inline  interval& operator=(interval I)
            {
                if (I._I!=this->_I)
                    mpri_swap(this->_I,I._I);
                return *this;
            }
            // inline  interval& operator=(const QQ& q)
            // {
            //     MPRI_SET_Q(_I,q.get_mpq_t());
            //     return *this;
            // }
        
            inline void set_l(const QQ&q)
            {
                mpq_set(mpri_lepref(this->_I),q.get_mpq_t());
            }
            inline void set_r(const QQ&q)
            {
                mpq_set(mpri_repref(this->_I),q.get_mpq_t());
            }

            inline QQ get_l() const
            {
                QQ ans;
                mpri_get_left(ans.get_mpq_t(),this->_I);
                return ans;
            }
            inline QQ get_r() const
            {
                QQ ans;
                mpri_get_right(ans.get_mpq_t(),this->_I);                
                return ans;
            }


            inline bool is_zero() const
            {
                return mpri_is_zero(_I);
            }
            inline bool is_notzero() const
            {
                return mpri_is_nonzero(_I);
            }
            inline bool is_nan() const
            {
                return mpria_mpq_is_nan(mpri_lepref(_I)) ||mpria_mpq_is_nan(mpri_repref(_I)); 
            }
            
            inline operator bool() const
            {
                return (this->is_notzero());
            }

            inline  bool operator==(const interval& I) const
            {
                // std::cout<<I<<std::endl;
                return mpri_equal(this->_I,I._I);
            }
            inline  bool operator==(int i) const
            {
                // std::cout<<I<<std::endl;
                // return mpri_equal(this->_I,I._I);
                return (mpq_cmp_si(mpri_lepref(this->_I),i,1)==0) &&
                       (mpq_cmp_si(mpri_repref(this->_I),i,1)==0) ;

            }
            inline  bool operator!=(int i) const
            {
                return !(*this==i);
            }
            inline  bool operator!=(const interval& I) const
            {
                return !(*this==I);
            }
            inline  interval operator+(const interval& I) const
            {
                interval ans;
                mpri_add(ans._I,this->_I,I._I);
                return ans;
            }
            inline  interval& operator+=(const interval& I) 
            {
                interval ans;
                mpri_add(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }
            
            
            inline  interval operator-(const interval& I) const
            {
                interval ans;
                mpri_sub(ans._I,this->_I,I._I);
                return ans;
            }
            inline interval operator&&(const interval& I) const
            {
                if (this->is_nan() || I.is_nan())
                    return interval::NAN_interval();
                // std::cout<<*this<<I<<std::endl;
                interval ans;
                {

                    if (QQ_cmp(mpri_lepref(this->_I),mpri_lepref(I._I))>=0)
                    {
                        mpq_set(mpri_lepref(ans._I),mpri_lepref(this->_I));
                    }
                    else{
                        mpq_set(mpri_lepref(ans._I),mpri_lepref(I._I));    
                    }
                    if (QQ_cmp(mpri_repref(this->_I),mpri_repref(I._I))<=0)
                    {
                        mpq_set(mpri_repref(ans._I),mpri_repref(this->_I));
                    }
                    else{
                        mpq_set(mpri_repref(ans._I),mpri_repref(I._I));    
                    }
                    if (QQ_cmp(mpri_lepref(ans._I),mpri_repref(ans._I))>0)
                    {
                        return interval::NAN_interval();
                    }

                }
                return ans; 
            }
            inline interval operator||(const interval& I) const
            {
                if (this->is_nan())
                    return I;
                if (I.is_nan())
                    return *this;
                interval ans;
                {

                    if (QQ_cmp(mpri_lepref(this->_I),mpri_lepref(I._I))<=0)
                    {
                        mpq_set(mpri_lepref(ans._I),mpri_lepref(this->_I));
                    }
                    else{
                        mpq_set(mpri_lepref(ans._I),mpri_lepref(I._I));    
                    }
                    if (QQ_cmp(mpri_repref(this->_I),mpri_repref(I._I))>=0)
                    {
                        mpq_set(mpri_repref(ans._I),mpri_repref(this->_I));
                    }
                    else{
                        mpq_set(mpri_repref(ans._I),mpri_repref(I._I));    
                    }
                    

                }
                return ans;
                // inline int QQ_cmp
                // return interval(std::min(this->get_l(),I.get_l()),std::max(this->get_r(),I.get_r()));           
            }
            inline static int QQ_cmp(const QQ& q1,const QQ& q2)// 0:q1=q2 1:q1>q2 -1:q1<q2
            {
                // int q1i=is_inf(q1);
                // int q2i=is_inf(q2);
                // if (q1i && q2i)
                // {
                //     if (q1i==q2i)
                //         return 0;
                //     if (q1i>q2i)
                //         return 1;
                //     return -1;
                // }
                // if (q1i)
                //     return q1i;
                // if (q2i)
                //     return -q2i;
                // return mpq_cmp(q1.get_mpq_t(),q2.get_mpq_t());
                return QQ_cmp(q1.get_mpq_t(),q2.get_mpq_t());
            }
            inline static int QQ_cmp(mpq_srcptr q1,mpq_srcptr q2)// 0:q1=q2 1:q1>q2 -1:q1<q2
            {
                int q1i=mpria_mpq_is_infinite(q1);
                int q2i=mpria_mpq_is_infinite(q2);
                if (q1i && q2i)
                {
                    if (q1i==q2i)
                        return 0;
                    if (q1i>q2i)
                        return 1;
                    return -1;
                }
                if (q1i)
                    return q1i;
                if (q2i)
                    return -q2i;
                return mpq_cmp(q1,q2);
            }
            inline  static int is_inf(const QQ & q)
            {
                return mpria_mpq_is_infinite(q.get_mpq_t());
            }
            inline interval operator-() const
            {
                interval ans;
                mpri_neg(ans._I,this->_I);
                return ans;
            }

            inline  interval& operator-=(const interval& I) 
            {
                interval ans;
                mpri_sub(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }


            inline  interval operator*(const interval& I) const
            {
                interval ans;
                mpri_mul(ans._I,this->_I,I._I);
                return ans;
            }
            inline  interval& operator*=(const interval& I) 
            {
                interval ans;
                mpri_mul(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }
            
            inline  interval operator/(const interval& I) const
            {
                interval ans;
                mpri_div(ans._I,this->_I,I._I);
                return ans;
            }
            inline  interval& operator/=(const interval& I) 
            {
                interval ans;
                mpri_div(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }
            
            
            inline interval sqr() const
            {
                interval ans;
                mpri_sqr(ans._I,this->_I);
                return ans;
            }
            
            inline interval pow(size_t n) const
            {
                if (n==1)
                    return *this;
                if (n==2)
                    return this->sqr();
                interval ans(1);
                if (n==0)
                    return ans;
                interval tmp(*this);
                size_t index=1;
                while (1){
                    if (index&n)
                    {
                        ans=ans*tmp;
                    }
                    index<<=1;
                    if (index>n)
                        break;
                    tmp=tmp.sqr();
                    if (index==n)
                        return tmp;
                }
                return ans;
            }
            
            void swap(interval & I)
            {   
                if (I._I!=this->_I)
                    mpri_swap(this->_I,I._I);
            }

            friend std::ostream& operator<<(std::ostream& stream, const interval& I)
            {
                _to_stream(_to_stream(stream<<"[",mpri_lepref(I._I))<<",",mpri_repref(I._I))<<"]";
                return stream;
            }

    

            inline static interval  inf_interval()
            {
                interval i;
                MPRIA_MPQ_SET_NEG_INF(mpri_lepref(i._I));
                MPRIA_MPQ_SET_POS_INF(mpri_repref(i._I));
                return i;
            } 
            inline static interval  NAN_interval()
            {
                interval i;
                MPRI_SET_NAN(i._I)
                return i;
            }
            inline static interval QQ_inf()
            {
                QQ q;
                MPRIA_MPQ_SET_POS_INF(q.get_mpq_t());
                return q;
            }
            inline static interval QQ_ninf()
            {
                QQ q;
                MPRIA_MPQ_SET_NEG_INF(q.get_mpq_t());
                return q;
            }
            

            
            
    };

    template <class comp=grlex>
    using interval_poly=polynomial_<interval,comp>;

    using interval_upoly=upolynomial_<interval>;

    template <class Tc, class To>
    interval_poly<To> assign(const polynomial_<Tc,To>& p,const std::map<variable,interval> & ass_list)
    {
        interval_poly<To> Pout(p.comp_ptr());
        basic_monomial<To> m(p.comp_ptr());
        basic_monomial<To> m1(p.comp_ptr());
        bool f=true;
        interval z, z1;
        for (auto &i:p)
        {

            m.clear();m.reserve(i.first.size());
            z=i.second;
            for (auto& j:i.first)
            {
                auto ptr=ass_list.find(j.first);
                if (ptr!=ass_list.end())
                {
                    z*=ptr->second.pow(j.second);
                }
                else
                    m.push_back(j);
            }
            if (!f && m!=m1)
            {
                if (z1)
                    Pout.push_back({std::move(m1),std::move(z1)});
                m1=std::move(m);
                z1=std::move(z);
            }
            else if (f)
            {
                m1=std::move(m);
                z1=std::move(z);
                f=false;
            }
            else
            {
                z1+=z;
            }

            
        }
        if (z1)
            Pout.push_back({std::move(m1),std::move(z1)});
        // std::cout<<Pout<<std::endl;
        Pout.normalization();
        return Pout;
    }



    template<class op_c>
    interval good_range(const upolynomial_<QQ> &p,const op_c& op) // op > || < >= <=
    {
        // std::cout<<"good_range:"<<p<<op.str()<<std::endl;
        const QQ zero=0;
        upolynomial_<ZZ> tmp_p;
        poly_convert(p,tmp_p);
        // std::cout<<"tmp_p:"<<tmp_p<<std::endl;
        upolynomial_<ZZ> G=tmp_p/polynomial_GCD(tmp_p,derivative(tmp_p));
        
        if (G.empty() ||  G.front().first.empty() && !op(G.front().second))
        {
            return interval::NAN_interval();
        }       
        if (G.front().first.empty())
        {
            interval I=interval::inf_interval();
            I.set_l(0);
            return I;
        }
        std::vector<std::pair<QQ,QQ>> l;
        // std::cout<<B<<std::endl;
        
        if (G.degree()==1)
        {
            
            // std::cout<<"G.degree()==1"<<std::endl;
            if (G.size()!=1)
            {
                QQ tmp=G.back().second;
                tmp=-tmp/G.front().second;
                // std::cout<<"tmp:"<<tmp<<std::endl;
                if (tmp>0)
                    l.push_back({tmp,tmp});
            }
        }
        else
        {
            ZZ B=RealRootBound(G);
            subuspensky(_upolynomial_Bto1(G,B),l,0,B);
        }
        if (l.empty())
        {
            if (!op(assign<QQ,ZZ,QQ>(p,1)))
            {
                return interval::NAN_interval();                
            }
            interval I=interval::inf_interval();
            I.set_l(0);
            return I;    
        }
        interval I=interval::inf_interval();
        I.set_l(0);
        QQ tmp_al,tmp_ar;
        tmp_al=assign<QQ,ZZ,QQ>(p,l.back().second);
        // std::cout<< l.back().second+1<<" "<< assign<QQ,ZZ,QQ>(p,l.back().second+1)<<" "<< op(assign<QQ,ZZ,QQ>(p,l.back().second+1))<<std::endl;
        
        if (!( op(assign<QQ,ZZ,QQ>(p,l.back().second+1))))
        {
            bool tmp=true;
            auto ptr=l.rbegin();
            for(;ptr!=l.rend();)
            {
                tmp_ar=std::move(tmp_al);
                tmp_al=assign<QQ,ZZ,QQ>(p,ptr->first);
                if (op(tmp_al) || op(zero) || op(tmp_ar)) 
                {
                    tmp=false;
                    I.set_r(ptr->second);
                    break;
                }
                auto ptr1=ptr;++ptr1;
                tmp_ar=std::move(tmp_al);
                if (ptr1!=l.rend())
                    tmp_al=assign<QQ,ZZ,QQ>(p,ptr->second);
                else 
                    tmp_al=assign<QQ,ZZ,QQ>(p,0);
                if (ptr1!=l.rend() && ( op(assign<QQ,ZZ,QQ>(p,(ptr1->second+ptr->first)/2))))
                {
                    tmp=false;
                    I.set_r(ptr->first);    
                    break;
                }
                ptr=ptr1;
            }
            if (tmp)
            {
                if (op(tmp_al)){
                    if (tmp_al==0 && !op(assign<QQ,ZZ,QQ>(p,l.front().first/2)))
                    {
                        I.set_r(0);
                        return I;    
                    }
                    I.set_r(l.front().first);
                    return I;
                }
                else
                    return  interval::NAN_interval();
            }
        }
        tmp_ar=assign<QQ,ZZ,QQ>(p,l.front().first);
        if (op(assign<QQ,ZZ,QQ>(p,0)) || op(assign<QQ,ZZ,QQ>(p,l.front().first/2)))
            I.set_l(0);
        else{
            bool tmp=true;
            
                for(auto ptr=l.begin();ptr!=l.end();)
                {
                    tmp_al=std::move(tmp_ar);
                    tmp_ar=assign<QQ,ZZ,QQ>(p,ptr->second);
                    if (op(tmp_al)|| op(tmp_ar) || op(zero))
                    {
                        tmp=false;
                    
                        I.set_l(ptr->first);
                        break;
                    }
                    auto ptr1=ptr;++ptr1;
                    tmp_al=std::move(tmp_ar);
                    if (ptr1!=l.end())
                        tmp_ar=assign<QQ,ZZ,QQ>(p,ptr1->first);
                    if (ptr1!=l.end() && op(assign<QQ,ZZ,QQ>(p,(ptr1->first+ptr->second)/2)))
                    {
                        tmp=false;
                    
                        I.set_l(ptr->second);    
                        break;
                    }
                    ptr=ptr1;
                }
                if (tmp)
                {
                    I.set_l(l.back().second);
                }
        }
        // std::cout<<I<<std::endl;
        return I;
    
    }
    interval feasible_range(const interval_upoly &p,char op);
    // {
    

    template <> 
    struct zore_check<interval>
    {
        inline bool operator()(const interval & I)
        {
            return I.is_zero();
        } 
    };

    template <class Tm,class To>
    std::ostream& operator<<  (std::ostream& stream, const basic_polynomial<Tm,interval,To>& v) 
    {
        if (v.size()==0)
            return stream<<'0';
        bool is_print=false;
        for(auto i=v.begin();i!=v.end();++i)
        {
            if (!zore_check<Tm>()(i->first) //&& !zore_check<Tc>()(i->second)
                ){
                if (is_print) //greater
                    stream<<"+";
                is_print=true;
                stream<< i->second<<"*";
                stream<<i->first;
                continue;
            }
            if (is_print) //greater
                stream<<"+";
            is_print=true;
            stream<<i->second;
        }
        return stream;
    }
}
#endif