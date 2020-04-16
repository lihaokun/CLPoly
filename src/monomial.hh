/*
Module Name:
    basic_monomial.hh
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
    template <class compare=lex>
    class basic_monomial
    {
        private:
            std::vector<std::pair<variable,int64_t>> __data; 
            int64_t __deg=0;       
        public:
            typedef compare compare_type;
            typedef variable first_type;
            typedef int64_t second_type;
            typedef typename std::vector<std::pair<variable,int64_t>>::iterator iterator;
            typedef typename std::vector<std::pair<variable,int64_t>>::const_iterator  const_iterator ;
            static const compare comp;
            constexpr int64_t deg() const {return this->__deg;}
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
                return tmp_deg==this->__deg && pair_vec_normal_check(this->__data,this->comp);
            }
            inline void normalization()
            {
                static std::function<bool(const std::pair<variable,int64_t> &,const std::pair<variable,int64_t> &)> tmp_comp=
                    [&](std::pair<variable,int64_t> a,std::pair<variable,int64_t> b){return (this->comp(a.first,b.first));};
                auto tmp_size=pair_vec_normalization(this->__data,tmp_comp);
                this->__data.resize(tmp_size);
                if (this->__data.size()*1.0/this->__data.capacity()<CLPOLY_shrink_bound && this->__data.capacity()>CLPOLY_shrink_bound_abs)
                {
                    //std::cout<<"shrink_to_fit";
                    this->__data.shrink_to_fit();
                }
                this->re_deg();
            }
            basic_monomial():__data(),__deg(0){}
            basic_monomial(const basic_monomial<compare> &m)
            :__data(m.__data),__deg(m.__deg)
            {}
            basic_monomial(const variable & v)
            :__data({{v,1}}),__deg(1)
            {}
            basic_monomial(basic_monomial<compare> &&m)
            :__data(std::move(m.__data)),__deg(std::move(m.__deg))
            {
                m.clear();
            }
            basic_monomial( std::initializer_list<std::pair<variable,int64_t>> init)
            :__data(init)
            {
                this->normalization();
            }
            basic_monomial(const std::vector<std::pair<variable,int64_t>> & v)
            :__data(v)
            {
                this->normalization();
            }
            basic_monomial(std::vector<std::pair<variable,int64_t>> && v)
            :__data(std::move(v))
            {
                this->normalization();
            }
            
            basic_monomial<compare> & operator=(const basic_monomial<compare> & m)
            {
                this->__data=m.__data;
                this->__deg=m.__deg;
                return *this;
            }
            basic_monomial<compare> & operator=(basic_monomial<compare> && m)
            {
                this->__data=std::move(m.__data);
                this->__deg=std::move(m.__deg);
                m.clear();
                return *this;
            }
            basic_monomial<compare> & operator=( std::initializer_list<std::pair<variable,int64_t>> init)
            {
                this->__data=init;
                this->normalization();
                return *this;
            }
            basic_monomial<compare> & operator=(const std::vector<std::pair<variable,int64_t>> & init)
            {
                this->__data=init;
                this->normalization();
                return *this;
            }
            basic_monomial<compare> & operator=(std::vector<std::pair<variable,int64_t>> && init)
            {
                this->__data=std::move(init);
                this->normalization();
                return *this;
            }

            constexpr auto size() const {return this->__data.size();}
            constexpr iterator begin() {return this->__data.begin();}
            constexpr const_iterator begin() const {return this->__data.begin();}
            
            constexpr iterator end() {return this->__data.end();}
            constexpr const_iterator end() const {return this->__data.end();}
            
            constexpr std::pair<variable,int64_t>& front() {return this->__data.front();}
            constexpr const std::pair<variable,int64_t>& front() const {return this->__data.front();}            
            constexpr std::pair<variable,int64_t>& back() {return this->__data.back();}
            constexpr const std::pair<variable,int64_t>& back() const {return this->__data.back();}
            constexpr void resize(std::size_t size){this->__data.resize(size);}
            constexpr void reserve(std::size_t size){this->__data.reserve(size);}
            inline void pop_back(){this->__deg-=this.back().second;this->__data.pop_back();}
            inline void push_back(const std::pair<variable,int64_t>&  value ){
                this->__deg+=value.second;
                this->__data.push_back(value);
            }
            inline void push_back(std::pair<variable,int64_t>&&  value){
                this->__deg+=value.second;
                this->__data.push_back(std::move(value));
            }
            constexpr std::pair<variable,int64_t>& operator[](std::size_t pos){return this->__data[pos];}
            constexpr const std::pair<variable,int64_t>& operator[](std::size_t pos) const {return this->__data[pos];}
            constexpr std::pair<variable,int64_t>& at(std::size_t pos){return this->__data.at(pos);}
            constexpr const std::pair<variable,int64_t>& at(std::size_t pos) const {return this->__data.at(pos);}
            constexpr bool empty()const{return this->__data.empty();} 
            constexpr std::vector<std::pair<variable,int64_t>> & data() {return this->__data;}
            constexpr const std::vector<std::pair<variable,int64_t>> & data() const {return this->__data;}
            inline void clear() 
            {
                this->__data.clear();
                this->__deg=0;
            }
            inline void swap(basic_monomial<compare> &m)
            {
                this->__data.swap(m.__data);
                std::swap(this->__deg,m.__deg);
            }

            inline basic_monomial<compare> operator*  (const basic_monomial<compare> &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left basic_monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right basic_monomial is not normal.(In left comparation.)");
                #endif
                basic_monomial<compare> new_p;
                new_p.__data=pair_vec_add(this->__data,p.__data,this->comp);
                new_p.__deg=this->__deg+p.__deg;
                return new_p;
            }
            inline basic_monomial<compare> operator/  (const basic_monomial<compare> &p)const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left basic_monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right basic_monomial is not normal.(In left comparation.)");
                #endif
                basic_monomial<compare> new_p;
                new_p.__data=pair_vec_sub(this->__data,p.__data,this->comp);
                new_p.__deg=this->__deg-p.__deg;
                return new_p;
            } 
            inline bool operator> (const basic_monomial<compare> &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left basic_monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right basic_monomial is not normal.(In left comparation.)");
                #endif
                if (this->deg()==p.deg())
                    return pair_vec_comp(this->__data,p.__data,this->comp);
                return this->deg()>p.deg();
            }
            inline bool operator== (const basic_monomial<compare> &p) const
            {
                #ifdef DEBUG
                    if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                        throw std::invalid_argument("Left basic_monomial is not normal.");
                    if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                        throw std::invalid_argument("Right basic_monomial is not normal.(In left comparation.)");
                #endif
                if (this->deg()==p.deg())
                    return pair_vec_equal_to(p.__data,this->__data);
                return false;
            }
            constexpr bool operator< (const basic_monomial<compare> &p) const
            {
                return (p>*this);
            }
            constexpr bool operator<=(const basic_monomial<compare> &p) const
            {
                return !(*this>p);
            }
            constexpr bool operator>=(const basic_monomial<compare> &p) const
            {
                return !(p>*this);
            }
            constexpr bool operator!=(const basic_monomial<compare> &p) const
            {
                return !(*this==p);
            }

            
            // inline basic_monomial operator-   () const
            // {
            //     return *this;
            // }   
            // inline basic_monomial operator+ () const
            // {
            //     basic_monomial new_m;
            //     new_m.reserve(this->size());
            //     for (auto &&i:*this)
            //     {
            //         new_m.push_back({i.first,negate(i.second)});
            //     }
            //     new_m.comp=this->comp;
            //     return new_m;
            // }   
            inline const_iterator find(const variable & m) const
            {
                std::function<bool(const std::pair<variable,int16_t> &,const std::pair<variable,int16_t> &)> tmp_comp=
                    [&](std::pair<variable,int16_t> a,std::pair<variable,int16_t> b){return (this->comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<variable,int16_t>(m,0),tmp_comp);
            }
            inline iterator find(const variable & m)
            {
                std::function<bool(const std::pair<variable,int16_t> &,const std::pair<variable,int16_t> &)> tmp_comp=
                    [&](std::pair<variable,int16_t> a,std::pair<variable,int16_t> b){return (this->comp(a.first,b.first));};
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
            friend std::ostream& operator<<  (std::ostream& stream, const basic_monomial<compare>& v) {
                if (v.empty())
                    return stream<<'1';
                bool is_print=false;
                for (const auto &i:v.data())
                { 
                    if (i.second!=0)
                    {
                        if (is_print)
                            stream<<"*";
                        stream<<i.first;
                        is_print=true;
                        if (i.second!=1)
                            stream<<"^"<<i.second;
                    }
                }
                return stream;
            }
            std::string  str() const{
                if (this->empty())
                    return "1";
                std::ostringstream ss;
                // bool is_print=false;
                // for (auto&&i:this->__data)
                // { 
                //     if (i.second!=0)
                //     {
                //         if (is_print)
                //             ss<<"*";
                //         ss<<i.first;
                //         is_print=true;
                //         if (i.second!=1)
                //             ss<<"^"<<i.second;
                //     }
                // }
                ss<<(*this);
                return ss.str();
            }

    };
    template<class compare> const compare basic_monomial<compare>::comp=compare();


    template<class compare>
    struct zore_check<basic_monomial<compare>>: public std::unary_function<basic_monomial<compare>, bool>
    {
        constexpr bool operator()(const basic_monomial<compare> & m)
        {
            return m.empty();
        } 
    };


    using monomial=basic_monomial<>;
    using grlex=graded<basic_monomial<lex>>;

}
namespace std{
    template<class compare>
    struct hash<clpoly::basic_monomial<compare>>: public std::unary_function<clpoly::basic_monomial<compare>, std::size_t>
    {
        std::size_t operator()(const clpoly::basic_monomial<compare> & m) const
        {
            return (hash<std::string>()(m.str()));
        }
    };
}
#endif