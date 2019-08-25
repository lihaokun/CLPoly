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
            std::vector<std::pair<variable,Tc>> __data;        
        public:
            static const std::function<bool(const variable &,const variable &)> init_comp;
            monomial():__data(){}
            monomial(const monomial<variable,Tc> &m)
            :__data(m.__data),__comp(m.__comp)
            {}
            monomial(monomial<variable,Tc> &&p)
            :__data(std::move(p.__data)),__comp(std::move(p.__comp))
            {
                p.__comp=init_comp;
            }
            atomic_polynomial( std::initializer_list<std::pair<variable,Tc>> init)
            :__data(init)
            {}
            atomic_polynomial(const std::vector<std::pair<variable,Tc>> & v)
            :__data(v)
            {}
            atomic_polynomial(std::vector<std::pair<variable,Tc>> && v)
            :__data(std::move(v))
            {}
            
            atomic_polynomial & operator=(const atomic_polynomial & p)
            {
                this->__comp=p.comp;
                this->__data=p.__data;
                return *this;
            }
            atomic_polynomial & operator=(atomic_polynomial && p)
            {
                this->__comp=std::move(p.comp);
                this->__data=std::move(p.__data);
                p.__comp=greater<variable>;
                return *this;
            }
            atomic_polynomial & operator=( std::initializer_list<std::pair<variable,Tc>> init)
            {
                this->__data=init;
                return *this;
            }
            inline auto size() const {return this->__data.size();}
            inline auto begin() {return this->__data.begin();}
            inline auto begin() const {return this->__data.begin();}
            
            inline auto end() {return this->__data.end();}
            inline auto end() const {return this->__data.end();}
            
            inline std::pair<variable,Tc>& back() {return this->__data.back();}
            inline const std::pair<variable,Tc>& back() const {return this->__data.back();}
            inline void resize(std::size_t size){this->__data.resize(size);}
            inline void reserve(std::size_t size){this->__data.reserve(size);}
            inline void pop_back(){this->__data.pop_back();}
            inline void push_back(const std::pair<variable,Tc>&  value ){this->__data.push_back(value);}
            inline void push_back(std::pair<variable,Tc>&&  value){this->__data.push_back(std::move(value));}
            inline std::pair<variable,Tc>& operator[](std::size_t pos){return this->__data[pos];}
            inline const std::pair<variable,Tc>& operator[](std::size_t pos) const {return this->__data.at(pos);}
            inline std::pair<variable,Tc>& at(std::size_t pos){return this->__data[pos];}
            inline const std::pair<variable,Tc>& at(std::size_t pos) const {return this->__data.at(pos);}
            inline bool empty()const{return this->__data.empty();} 
            inline void clear() 
            {
                this->__data.clear();
                this->__comp=greater<variable>;
            }
            inline void swap(atomic_polynomial<variable,Tc> &p)
            {
                this->__data.swap(p.__data);
                this->__comp.swap(p.__comp);
            }

            
            
            inline bool is_normal()  const 
            {
                return pair_vec_normal_check(this->begin(),this->end(),this->__comp);
            }
            inline void normalization()
            {
                std::function<bool(const std::pair<variable,Tc> &,const std::pair<variable,Tc> &)> tmp_comp=
                    [&](std::pair<variable,Tc> a,std::pair<variable,Tc> b){return (this->__comp(a.first,b.first));};
                auto tmp_size=pair_vec_normalization(this->begin(),this->end(),tmp_comp);
                this->resize(tmp_size);
            }
            inline atomic_polynomial<variable,Tc> operator+  (const atomic_polynomial<variable,Tc> &p)const
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
            inline atomic_polynomial<variable,Tc> operator-  (const atomic_polynomial<variable,Tc> &p)const
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
            inline atomic_polynomial<variable,Tc> operator+  () const
            {
                return *this;
            }   
            inline atomic_polynomial<variable,Tc> operator- () const
            {
                atomic_polynomial<variable,Tc> new_p;
                new_p.reserve(this->size());
                for (auto &i:*this)
                {
                    new_p.push_back({i.first,negate(i.second)});
                }
                new_p.__comp=this->__comp;
                return new_p;
            }   
    };
}
#endif