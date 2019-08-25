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
namespace clpoly{
    class monomial
    {
        private:
            std::function<bool(const variable &,const variable &)> __comp=init_comp;
            std::vector<std::pair<variable,int64_t>> __data; 
            int64_t __deg=0;       
        public:
            static const std::function<bool(const variable &,const variable &)> init_comp;
            inline int64_t deg(){return this->__deg;}
            inline int64_t re_deg()
            {
                this->__deg=0;
                for(auto &i:this->__data)
                    this->__deg+=i.second;
                return this->__deg;
            }
            inline bool is_normal()  const 
            {
                int64_t tmp_deg=0;
                for(auto &i:this->__data)
                    tmp_deg+=i.second;
                return tmp_deg==this->__deg && pair_vec_normal_check(this->__data.begin(),this->__data.end(),this->__comp);
            }
            inline void normalization()
            {
                std::function<bool(const std::pair<variable,int64_t> &,const std::pair<variable,int64_t> &)> tmp_comp=
                    [&](std::pair<variable,int64_t> a,std::pair<variable,int64_t> b){return (this->__comp(a.first,b.first));};
                auto tmp_size=pair_vec_normalization(this->__data.begin(),this->__data.end(),tmp_comp);
                this->__data.resize(tmp_size);
                this->re_deg();
            }
            monomial():__data(),__deg(0){}
            monomial(const monomial &m)
            :__data(m.__data),__comp(m.__comp),__deg(m.__deg)
            {}
            monomial(monomial &&m)
            :__data(std::move(m.__data)),__comp(std::move(m.__comp)),__deg(std::move(m.__deg))
            {
                m.__comp=init_comp;
                m.__deg=0;

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
                m.__comp=init_comp;
                return *this;
            }
            monomial & operator=( std::initializer_list<std::pair<variable,int64_t>> init)
            {
                this->__data=init;
                this->normalization();
                return *this;
            }
            inline auto size() const {return this->__data.size();}
            inline auto begin() {return this->__data.begin();}
            inline auto begin() const {return this->__data.begin();}
            
            inline auto end() {return this->__data.end();}
            inline auto end() const {return this->__data.end();}
            
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
            inline const std::pair<variable,int64_t>& operator[](std::size_t pos) const {return this->__data.at(pos);}
            inline std::pair<variable,int64_t>& at(std::size_t pos){return this->__data[pos];}
            inline const std::pair<variable,int64_t>& at(std::size_t pos) const {return this->__data.at(pos);}
            inline bool empty()const{return this->__data.empty();} 
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
            inline std::function<bool(const variable &,const variable &)> & comp()  {return this->__comp;}
            
            
            inline monomial operator+  (const monomial &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                auto new_p=pair_vec_add(*this,p,this->__comp);
                new_p.__comp=this->__comp;
                return new_p;
            }
            inline monomial operator-  (const monomial &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->__comp))
                        throw std::invalid_argument("Left atomic_polynomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->__comp))
                        throw std::invalid_argument("Right atomic_polynomial is not normal.(In left comparation.)");
                #endif
                auto new_p=pair_vec_sub(*this,p,this->__comp);
                new_p.__comp=this->__comp;
                return new_p;
            } 
            inline monomial operator+  () const
            {
                return *this;
            }   
            inline monomial operator- () const
            {
                monomial new_m;
                new_m.reserve(this->size());
                for (auto &i:*this)
                {
                    new_m.push_back({i.first,negate(i.second)});
                }
                new_m.__comp=this->__comp;
                return new_m;
            }   
    };
}
#endif