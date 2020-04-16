/*
Module Name:
    atomic_polynomial.hh
Abstract:
    定义类：atomic_polynomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_ATOMIC_POLYNOMIAL_HH
#define CLPOLY_ATOMIC_POLYNOMIAL_HH
#include "basic.hh"
#include <vector>
#include <functional>
#include <sstream>
namespace clpoly{

    template <class Tm,class Tc,class compare=greater<Tm>>
    class atomic_polynomial
    {
        private:
            std::vector<std::pair<Tm,Tc>> __data;
            
        public:
            typedef compare compare_type;
            typedef Tm momomial_type;
            typedef Tc coeff_type;
            typedef Tm first_type;
            typedef Tc second_type;

            typedef typename std::vector<std::pair<Tm,Tc>>::iterator iterator;
            typedef typename std::vector<std::pair<Tm,Tc>>::const_iterator  const_iterator ;
            
            static const compare comp;
            static const Tc Tc_zero;
            atomic_polynomial():__data(){}
            atomic_polynomial(const atomic_polynomial<Tm,Tc,compare> &p)
            :__data(p.__data)
            {}
            atomic_polynomial(atomic_polynomial<Tm,Tc> &&p)
            :__data(std::move(p.__data))
            {}
            atomic_polynomial( std::initializer_list<std::pair<Tm,Tc>> init)
            :__data(init)
            {
                this->normalization();
            }
            atomic_polynomial(const std::vector<std::pair<Tm,Tc>> & v)
            :__data(v)
            {
                this->normalization();
            }
            atomic_polynomial(std::vector<std::pair<Tm,Tc>> && v)
            :__data(std::move(v))
            {
                this->normalization();
            }
            
            inline atomic_polynomial<Tm,Tc,compare> & operator=(const atomic_polynomial<Tm,Tc> & p)
            {
                this->__data=p.__data;
                return *this;
            }
            
            inline atomic_polynomial<Tm,Tc,compare> & operator=(atomic_polynomial<Tm,Tc> && p)
            {
                this->__data=std::move(p.__data);
                return *this;
            }

            inline atomic_polynomial<Tm,Tc,compare> & operator=( std::initializer_list<std::pair<Tm,Tc>> init)
            {
                this->__data=init;
                this->normalization();
                return *this;
            }
            
            inline atomic_polynomial<Tm,Tc,compare> & operator=(const std::vector<std::pair<Tm,Tc>> & init)
            {
                this->__data=init;
                this->normalization();
                return *this;
            }

            inline atomic_polynomial<Tm,Tc,compare> & operator=(std::vector<std::pair<Tm,Tc>> && init)
            {
                this->__data=std::move(init);
                this->normalization();
                return *this;
            }

            constexpr auto size() const {return this->__data.size();}
            constexpr iterator begin() {return this->__data.begin();}
            constexpr const_iterator  begin() const {return this->__data.begin();}
            constexpr iterator  end() {return this->__data.end();}
            constexpr const_iterator  end() const {return this->__data.end();}
            constexpr std::pair<Tm,Tc>& front() {return this->__data.front();}
            constexpr const std::pair<Tm,Tc>& front() const {return this->__data.front();}
            constexpr std::pair<Tm,Tc>& back() {return this->__data.back();}
            constexpr const std::pair<Tm,Tc>& back() const {return this->__data.back();}
            constexpr void resize(std::size_t size){this->__data.resize(size);}
            constexpr void reserve(std::size_t size){this->__data.reserve(size);}
            constexpr void pop_back(){this->__data.pop_back();}
            constexpr void push_back(const std::pair<Tm,Tc>&  value ){this->__data.push_back(value);}
            constexpr void push_back(std::pair<Tm,Tc>&&  value){this->__data.push_back(std::move(value));}
            constexpr std::pair<Tm,Tc>& operator[](std::size_t pos){return this->__data[pos];}
            constexpr const std::pair<Tm,Tc>& operator[](std::size_t pos) const {return this->__data[pos];}
            constexpr std::pair<Tm,Tc>& at(std::size_t pos){return this->__data.at(pos);}
            constexpr const std::pair<Tm,Tc>& at(std::size_t pos) const {return this->__data.at(pos);}
            constexpr bool empty()const{return this->__data.empty();} 
            constexpr std::vector<std::pair<Tm,Tc>> & data() {return this->__data;}
            constexpr const std::vector<std::pair<Tm,Tc>> & data() const {return this->__data;}
            
            constexpr void shrink_to_fit(){this->__data.shrink_to_fit();}
            constexpr std::size_t capacity(){return this->__data.capacity();}
            inline void clear() 
            {
                this->__data.clear();
            }

            inline void swap(atomic_polynomial<Tm,Tc> &p)
            {
                this->__data.swap(p.__data);
            }

            
            
            inline bool is_normal()  const 
            {
                return pair_vec_normal_check(this->__data,this->comp);
                // if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //     return false
                // for (auto && i:*this)
                //     if (!i.first.is_normal())
                //         return false;
                // return true;
            
            }

            inline void normalization()
            {
                // for (auto && i:*this)
                //     i.first.normalization());
                static std::function<bool(const std::pair<Tm,Tc> &,const std::pair<Tm,Tc> &)> tmp_comp=
                    [&](std::pair<Tm,Tc> a,std::pair<Tm,Tc> b){return (this->comp(a.first,b.first));};
                auto tmp_size=pair_vec_normalization(this->__data,tmp_comp);
                this->resize(tmp_size);
                if (this->__data.size()*1.0/this->__data.capacity()<CLPOLY_shrink_bound && this->__data.size()>CLPOLY_shrink_bound_abs)
                {
                    //std::cout<<"shrink_to_fit";
                    this->__data.shrink_to_fit();
                }

            }

            inline atomic_polynomial<Tm,Tc> operator+  (const atomic_polynomial<Tm,Tc> &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                auto new_p=pair_vec_add(this->__data,p.__data,this->comp);
                return new_p;
            }

            inline atomic_polynomial<Tm,Tc> operator-  (const atomic_polynomial<Tm,Tc> &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                auto new_p=pair_vec_sub(this->__data,p.__data,this->comp);
                return new_p;
            } 
            
            inline atomic_polynomial<Tm,Tc> operator+  () const
            {
                return *this;
            }   

            inline atomic_polynomial<Tm,Tc> operator- () const
            {
                atomic_polynomial<Tm,Tc> new_p;
                new_p.__data=pair_vec_negate(this->__data);
                return new_p;
            } 

            inline atomic_polynomial<Tm,Tc> operator*(const atomic_polynomial<Tm,Tc> & p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                atomic_polynomial<Tm,Tc> new_p;
                new_p.__data= pair_vec_multiplies(this->__data,p.__data,this->comp);    
                return new_p;
            }
            inline atomic_polynomial<Tm,Tc>  power(unsigned i) const
            {
                atomic_polynomial<Tm,Tc> newp;
                if (this->empty())
                    return newp; 
                atomic_polynomial<Tm,Tc> p(*this);
                while (i!=0)
                {
                    if (i%2!=0)
                        if (newp.empty())
                            newp=p;
                        else
                            newp=newp*p;
                    i>>=1;
                    if (i!=0)
                        p=p*p;
                }
                return newp;
            }
            inline bool operator> (const atomic_polynomial<Tm,Tc> &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                return pair_vec_comp(this->__data,p.__data,this->comp);
            }
            inline bool operator< (const atomic_polynomial<Tm,Tc> &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                return pair_vec_comp(p.__data,this->__data,this->comp);
            }
            inline bool operator== (const atomic_polynomial<Tm,Tc> &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                return pair_vec_equal_to(p.__data,this->__data);
            }
            inline bool operator<=(const atomic_polynomial<Tm,Tc> &p) const
            {
                return !(*this>p);
            }
            inline bool operator>=(const atomic_polynomial<Tm,Tc> &p) const
            {
                return !(*this<p);
            }
            std::string str() const
            {
                if (this->size()==0)
                    return "0";
                std::ostringstream ss;
                // bool is_print=false;
                // for(auto i=this->begin();i!=this->end();++i)
                // {
                //     if (!zore_check(i->first) && !zore_check(i->second)){
                //         if (is_print && greater(i->second,Tc_zero))
                //             ss<<"+";
                //         is_print=true;
                //         if (equal_to(i->second,-1))
                //             ss<<"-";
                //         else
                //             if (!equal_to(i->second,1) )
                //                 ss<< i->second<<"*";
                //         ss<<i->first;
                //         continue;
                //     }
                //     if (zore_check(i->first) && !zore_check(i->second))
                //         if (is_print && greater(i->second,Tc_zero))
                //             ss<<"+";
                //         is_print=true;
                //         ss<<i->second;
                // }
                ss<<(*this);
                return ss.str();
            }
            inline const_iterator find(const Tm & m) const
            {
                std::function<bool(const std::pair<Tm,Tc> &,const std::pair<Tm,Tc> &)> tmp_comp=
                    [&](std::pair<Tm,Tc> a,std::pair<Tm,Tc> b){return (this->comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<Tm,Tc>(m,0),tmp_comp);
            }
            inline iterator find(const Tm & m)
            {
                std::function<bool(const std::pair<Tm,Tc> &,const std::pair<Tm,Tc> &)> tmp_comp=
                    [&](std::pair<Tm,Tc> a,std::pair<Tm,Tc> b){return (this->comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<Tm,Tc>(m,0),tmp_comp);
            }
            inline const Tc & coef(const Tm & m) const
            {
                auto l=this->find(m);
                if (l!=this->end())
                    return l->second;
                else
                {
                    //Tc tmp_Tc=0;
                    return Tc_zero;//tmp_Tc;
                }
                    
            }
            friend std::ostream& operator<<  (std::ostream& stream, const atomic_polynomial<Tm,Tc>& v) 
            {
                if (v.size()==0)
                    return stream<<'0';
                bool is_print=false;
                for(auto i=v.begin();i!=v.end();++i)
                {
                    if (!zore_check<Tm>()(i->first) && !zore_check<Tc>()(i->second)){
                        if (is_print && i->second>Tc_zero) //greater
                            stream<<"+";
                        is_print=true;
                        if (i->second==-1) //equal_to
                            stream<<"-";
                        else
                            if (i->second!=1) //equal_to
                                stream<< i->second<<"*";
                        stream<<i->first;
                        continue;
                    }
                    if (zore_check<Tm>()(i->first) && !zore_check<Tc>()(i->second))
                        if (is_print && i->second>Tc_zero) //greater
                            stream<<"+";
                        is_print=true;
                        stream<<i->second;
                }
                return stream;
            }
              
    };

    template<class Tm,class Tc,class compare> const compare  atomic_polynomial<Tm, Tc, compare>::comp=compare();
    template<class Tm,class Tc,class compare> const Tc  atomic_polynomial<Tm, Tc, compare>::Tc_zero=0;
   
}
namespace std{
    template<class Tm,class Tc>
    struct hash<clpoly::atomic_polynomial<Tm,Tc>>
    {
        std::size_t operator()(const clpoly::atomic_polynomial<Tm,Tc> & m) const
        {
            return (hash<std::string>()(m.str()));
        }
    };
}
#endif