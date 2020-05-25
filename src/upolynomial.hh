/*
Module Name:
    upolynomial.hh
Abstract:
    定义upolynomial
Author:
    haokun li
Notes:
    废弃
*/
#ifndef CLPOLY_UPOLYNOMIAL_HH
#define CLPOLY_UPOLYNOMIAL_HH
#include "polynomial.hh"
#include "basic_polynomial.hh"
#include "variable.hh"
#include <vector>
namespace clpoly{
    template <class coeff>
    class upolynomial
    {
    private:
        std::vector<std::pair<int64_t,coeff>> __data;
        variable __var;
    public:
        typedef coeff coeff_type;
        typedef typename std::vector<std::pair<int64_t,coeff>>::iterator iterator;
        typedef typename std::vector<std::pair<int64_t,coeff>>::const_iterator  const_iterator ;
            
        upolynomial()==default;
        upolynomial(const upolynomial<coeff> & p)==default;
        upolynomial(upolynomial<coeff> && p)==default;
        
        
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
        constexpr int64_t deg() const {return __data.front().first;}
        constexpr const variable & var() const {return __var;}

        void clear()
        {
            __data.clear();
        }
    };
    
    
    
    
}
#endif