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
namespace clpoly{

    template <class Tm,class Tc>
    class atomic_polynomial
    {
        private:
            std::function<bool(const Tm &,const Tm &)> __comp=init_comp;
            std::vector<std::pair<Tm,Tc>> __data;
        public:
            static const std::function<bool(const Tm &,const Tm &)> init_comp;
            atomic_polynomial():__data(){}
            atomic_polynomial(const atomic_polynomial<Tm,Tc> &p)
            :__data(p.__data),__comp(p.__comp)
            {}
            atomic_polynomial(atomic_polynomial<Tm,Tc> &&p)
            :__data(std::move(p.__data)),__comp(std::move(p.__comp))
            {
                p.__comp=init_comp;
            }
            atomic_polynomial( std::initializer_list<std::pair<Tm,Tc>> init)
            :__data(init)
            {}
            atomic_polynomial(const std::vector<std::pair<Tm,Tc>> & v)
            :__data(v)
            {}
            atomic_polynomial(std::vector<std::pair<Tm,Tc>> && v)
            :__data(std::move(v))
            {}
            
            atomic_polynomial<Tm,Tc> & operator=(const atomic_polynomial<Tm,Tc> & p)
            {
                this->__comp=p.__comp;
                this->__data=p.__data;
                return *this;
            }
            atomic_polynomial<Tm,Tc> & operator=(atomic_polynomial<Tm,Tc> && p)
            {
                this->__comp=std::move(p.__comp);
                this->__data=std::move(p.__data);
                p.__comp=init_comp;
                return *this;
            }
            atomic_polynomial<Tm,Tc> & operator=( std::initializer_list<std::pair<Tm,Tc>> init)
            {
                this->__data=init;
                return *this;
            }
            inline auto size() const {return this->__data.size();}
            inline auto begin() {return this->__data.begin();}
            inline auto begin() const {return this->__data.begin();}
            
            inline auto end() {return this->__data.end();}
            inline auto end() const {return this->__data.end();}
            
            inline std::pair<Tm,Tc>& back() {return this->__data.back();}
            inline const std::pair<Tm,Tc>& back() const {return this->__data.back();}
            inline void resize(std::size_t size){this->__data.resize(size);}
            inline void reserve(std::size_t size){this->__data.reserve(size);}
            inline void pop_back(){this->__data.pop_back();}
            inline void push_back(const std::pair<Tm,Tc>&  value ){this->__data.push_back(value);}
            inline void push_back(std::pair<Tm,Tc>&&  value){this->__data.push_back(std::move(value));}
            inline std::pair<Tm,Tc>& operator[](std::size_t pos){return this->__data[pos];}
            inline const std::pair<Tm,Tc>& operator[](std::size_t pos) const {return this->__data.at(pos);}
            inline std::pair<Tm,Tc>& at(std::size_t pos){return this->__data[pos];}
            inline const std::pair<Tm,Tc>& at(std::size_t pos) const {return this->__data.at(pos);}
            inline bool empty()const{return this->__data.empty();} 
            inline std::pair<Tm,Tc>* data() {return this->__data.data();}
            inline const std::pair<Tm,Tc>* data() const {return this->__data.data();}
            inline const std::function<bool(const Tm &,const Tm &)> & comp() const {return this->__comp;}
            inline std::function<bool(const Tm &,const Tm &)> & comp()  {return this->__comp;}
            

            inline void clear() 
            {
                this->__data.clear();
                this->__comp=init_comp;
            }
            inline void swap(atomic_polynomial<Tm,Tc> &p)
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
                std::function<bool(const std::pair<Tm,Tc> &,const std::pair<Tm,Tc> &)> tmp_comp=
                    [&](std::pair<Tm,Tc> a,std::pair<Tm,Tc> b){return (this->__comp(a.first,b.first));};
                auto tmp_size=pair_vec_normalization(this->begin(),this->end(),tmp_comp);
                this->resize(tmp_size);
            }
            inline atomic_polynomial<Tm,Tc> operator+  (const atomic_polynomial<Tm,Tc> &p)const
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
            inline atomic_polynomial<Tm,Tc> operator-  (const atomic_polynomial<Tm,Tc> &p)const
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
            inline atomic_polynomial<Tm,Tc> operator+  () const
            {
                return *this;
            }   
            inline atomic_polynomial<Tm,Tc> operator- () const
            {
                atomic_polynomial<Tm,Tc> new_p;
                new_p.reserve(this->size());
                for (auto &i:*this)
                {
                    new_p.push_back({i.first,negate(i.second)});
                }
                new_p.__comp=this->__comp;
                return new_p;
            }   
    };
    template<class Tm,class Tc> const std::function<bool(const Tm &,const Tm &)>  atomic_polynomial<Tm,Tc>::init_comp(greater<Tm>);
}
#endif