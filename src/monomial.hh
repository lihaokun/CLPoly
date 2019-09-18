/*
Module Name:
    monomial.hh
Abstract:
    定义类：monomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_MONOMIAL_HH
#define CLPOLY_MONOMIAL_HH
#include "variable.hh"
#include "basic.hh"
#include <vector>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
namespace clpoly{
    class monomial
    {
        private:
            std::function<bool(const variable &,const variable &)> __comp=init_comp;
            std::vector<std::pair<variable,int64_t>> __data; 
            int64_t __deg=0;       
        public:
            typedef typename std::vector<std::pair<variable,int64_t>>::iterator iterator;
            typedef typename std::vector<std::pair<variable,int64_t>>::const_iterator  const_iterator ;
            static const std::function<bool(const variable &,const variable &)> init_comp;
            inline int64_t deg() const {return this->__deg;}
            inline int64_t re_deg()
            {
                this->__deg=0;
                for(auto &&i:this->__data)
                    this->__deg+=i.second;
                return this->__deg;
            }
            inline bool is_normal()  const 
            {
                int64_t tmp_deg=0;
                for(auto &&i:this->__data)
                    tmp_deg+=i.second;
                return tmp_deg==this->__deg && pair_vec_normal_check(this->__data.begin(),this->__data.end(),this->__comp);
            }
            inline void normalization()
            {
                std::function<bool(const std::pair<variable,int64_t> &,const std::pair<variable,int64_t> &)> tmp_comp=
                    [&](std::pair<variable,int64_t> a,std::pair<variable,int64_t> b){return (this->__comp(a.first,b.first));};
                auto tmp_size=pair_vec_normalization(this->__data.begin(),this->__data.end(),tmp_comp);
                this->__data.resize(tmp_size);
                if (this->__data.size()*1.0/this->__data.capacity()<shrink_bound && this->__data.size()>shrink_bound_abs)
                {
                    //std::cout<<"shrink_to_fit";
                    this->__data.shrink_to_fit();
                }
                this->re_deg();
            }
            monomial():__data(),__deg(0){}
            monomial(const monomial &m)
            :__data(m.__data),__comp(m.__comp),__deg(m.__deg)
            {}
            monomial(const variable & v)
            :__data({{v,1}}),__deg(1)
            {}
            monomial(monomial &&m)
            :__data(std::move(m.__data)),__comp(std::move(m.__comp)),__deg(std::move(m.__deg))
            {
                // m.__comp=init_comp;
                // m.__deg=0;
                m.clear();

            }
            monomial( std::initializer_list<std::pair<variable,int64_t>> init)
            :__data(init)
            {
                this->normalization();
            }
            monomial(const std::vector<std::pair<variable,int64_t>> & v,const std::function<bool(const variable &,const variable &)>& comp=init_comp)
            :__data(v),__comp(comp)
            {
                this->normalization();
            }
            monomial(std::vector<std::pair<variable,int64_t>> && v,const std::function<bool(const variable &,const variable &)>& comp=init_comp)
            :__data(std::move(v)),__comp(comp)
            {
                this->normalization();
            }
            monomial(std::vector<std::pair<variable,int64_t>> && v,std::function<bool(const variable &,const variable &)>&& comp)
            :__data(std::move(v)),__comp(std::move(comp))
            {
                this->normalization();
            }
            
            
            monomial & operator=(const monomial & m)
            {
                this->__comp=m.__comp;
                this->__data=m.__data;
                this->__deg=m.__deg;
                return *this;
            }
            monomial & operator=(monomial && m)
            {
                this->__comp=std::move(m.__comp);
                this->__data=std::move(m.__data);
                this->__deg=std::move(m.__deg);
                m.clear();
                //m.__comp=init_comp;
                return *this;
            }
            monomial & operator=( std::initializer_list<std::pair<variable,int64_t>> init)
            {
                this->__data=init;
                this->normalization();
                return *this;
            }
            inline auto size() const {return this->__data.size();}
            inline iterator begin() {return this->__data.begin();}
            inline const_iterator begin() const {return this->__data.begin();}
            
            inline iterator end() {return this->__data.end();}
            inline const_iterator end() const {return this->__data.end();}
            
            inline std::pair<variable,int64_t>& back() {return this->__data.back();}
            inline const std::pair<variable,int64_t>& back() const {return this->__data.back();}
            inline void resize(std::size_t size){this->__data.resize(size);}
            inline void reserve(std::size_t size){this->__data.reserve(size);}
            inline void pop_back(){this->__data.pop_back();}
            inline void push_back(const std::pair<variable,int64_t>&  value ){
                this->__deg+=value.second;
                this->__data.push_back(value);
            }
            inline void push_back(std::pair<variable,int64_t>&&  value){
                this->__deg+=value.second;
                this->__data.push_back(std::move(value));
            }
            inline std::pair<variable,int64_t>& operator[](std::size_t pos){return this->__data[pos];}
            inline const std::pair<variable,int64_t>& operator[](std::size_t pos) const {return this->__data[pos];}
            inline std::pair<variable,int64_t>& at(std::size_t pos){return this->__data.at(pos);}
            inline const std::pair<variable,int64_t>& at(std::size_t pos) const {return this->__data.at(pos);}
            inline bool empty()const{return this->__data.empty();} 
            inline std::vector<std::pair<variable,int64_t>> & data() {return this->__data;}
            inline const std::vector<std::pair<variable,int64_t>> & data() const {return this->__data;}
            inline void clear() 
            {
                this->__data.clear();
                this->__comp=init_comp;
                this->__deg=0;
            }
            inline void swap(monomial &m)
            {
                this->__data.swap(m.__data);
                this->__comp.swap(m.__comp);
                std::swap(this->__deg,m.__deg);
            }
            inline const std::function<bool(const variable &,const variable &)> & comp() const {return this->__comp;}
            inline void comp(const std::function<bool(const variable &,const variable &)> & c)  {this->__comp=c;this->normalization();}
            inline void comp(std::function<bool(const variable &,const variable &)> && c)  {this->__comp=std::move(c);this->normalization();}
            
            
            inline monomial operator*  (const monomial &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right monomial is not normal.(In left comparation.)");
                #endif
                monomial new_p;
                new_p.__data=pair_vec_add(this->__data,p.__data,this->__comp);
                new_p.__comp=this->__comp;
                new_p.__deg=this->__deg+p.__deg;
                return new_p;
            }
            inline monomial operator/  (const monomial &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right monomial is not normal.(In left comparation.)");
                #endif
                monomial new_p;
                new_p.__data=pair_vec_sub(this->__data,p.__data,this->__comp);
                new_p.__comp=this->__comp;
                new_p.__deg=this->__deg-p.__deg;
                return new_p;
            } 
            // inline monomial operator*  (const monomial &p)const {return *this+p;}
            // inline monomial operator/  (const monomial &p)const {return *this-p;}
            inline bool operator> (const monomial &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right monomial is not normal.(In left comparation.)");
                #endif
                if (this->deg()==p.deg())
                    return pair_vec_comp(this->__data,p.__data,this->__comp);
                return this->deg()>p.deg();
            }
            inline bool operator< (const monomial &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right monomial is not normal.(In left comparation.)");
                #endif
                if (this->deg()==p.deg())
                    return pair_vec_comp(p.__data,this->__data,this->__comp);
                return this->deg()>p.deg();
            }
            inline bool operator== (const monomial &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right monomial is not normal.(In left comparation.)");
                #endif
                if (this->deg()==p.deg())
                    return pair_vec_equal_to(p.__data,this->__data);
                return false;
            }
            inline bool operator<=(const monomial &p) const
            {
                return !(*this>p);
            }
            inline bool operator>=(const monomial &p) const
            {
                return !(*this<p);
            }
            
            
            
            // inline monomial operator-   () const
            // {
            //     return *this;
            // }   
            // inline monomial operator+ () const
            // {
            //     monomial new_m;
            //     new_m.reserve(this->size());
            //     for (auto &&i:*this)
            //     {
            //         new_m.push_back({i.first,negate(i.second)});
            //     }
            //     new_m.__comp=this->__comp;
            //     return new_m;
            // }   
            inline const_iterator find(const variable & m) const
            {
                std::function<bool(const std::pair<variable,int16_t> &,const std::pair<variable,int16_t> &)> tmp_comp=
                    [&](std::pair<variable,int16_t> a,std::pair<variable,int16_t> b){return (this->__comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<variable,int16_t>(m,0),tmp_comp);
            }
            inline iterator find(const variable & m)
            {
                std::function<bool(const std::pair<variable,int16_t> &,const std::pair<variable,int16_t> &)> tmp_comp=
                    [&](std::pair<variable,int16_t> a,std::pair<variable,int16_t> b){return (this->__comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<variable,int16_t>(m,0),tmp_comp);
            }
            inline int64_t deg(const variable & m) const
            {
                auto l=this->find(m);
                if (l!=this->end())
                    return l->second;
                else
                    return 0;
            }
            std::string  str() const;
            friend std::ostream& operator<<  (std::ostream& stream, const monomial& v) {
                return stream<<v.str();
            }
    };
    template<>
    inline bool zore_check(const monomial& m)
    {
        return m.empty();
    }
}
namespace std{
    template<>
    struct hash<clpoly::monomial>
    {
        std::size_t operator()(const clpoly::monomial & m) const
        {
            return (hash<std::string>()(m.str()));
        }
    };
}
#endif