/*
Module Name:
    basic_polynomial.hh
Abstract:
    定义类：basic_polynomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_ATOMIC_POLYNOMIAL_HH
#define CLPOLY_ATOMIC_POLYNOMIAL_HH
#include <vector>
#include <list>
#include <functional>
#include <algorithm>
#include <sstream>
#include <clpoly/variable.hh>
namespace clpoly{
    class basic_polynomial_status{
        public:
            bool is_var;
            bool is_deg;
            std::list<std::pair<variable,std::int64_t>>  variables;
            int64_t deg;

            basic_polynomial_status()
            :is_var(false),is_deg(false)
            {}

            basic_polynomial_status(const basic_polynomial_status & b)
            :is_var(b.is_var),is_deg(b.is_deg),deg(b.deg)
            {
                if (is_var)
                    variables=b.variables;
            }

            basic_polynomial_status(basic_polynomial_status && b)
            :is_var(std::move(b.is_var)),is_deg(std::move(b.is_deg)),deg(std::move(b.deg))
            {
                if (is_var)
                    variables=std::move(b.variables);
                b.clear();
            }

            basic_polynomial_status & operator=(const basic_polynomial_status & b)
            {
                if (&b==this)
                    return *this;
                this->is_var=b.is_var;
                this->deg=b.deg;
                this->is_deg=b.is_deg;
                if (this->is_var)
                    this->variables=b.variables;
                return *this;
            }

            basic_polynomial_status & operator=(basic_polynomial_status && b)
            {
                if (&b==this)
                    return *this;
                this->is_var=b.is_var;
                this->deg=b.deg;
                this->is_deg=b.is_deg;
                if (this->is_var)
                    this->variables=std::move(b.variables);
                b.clear();
                return *this;
            }

            void swap(basic_polynomial_status & b)
            {
                if (&b==this)
                    return void();
                std::swap(this->is_deg,b.is_deg);
                std::swap(this->is_var,b.is_var);
                std::swap(this->deg,b.deg);
                if (this->is_var || b.is_var)
                    std::swap(this->variables,b.variables);
            }
            void clear()
            {
                is_var=is_deg=false;
            }
    };
    template <class Tm,class Tc,class compare>
    class basic_polynomial
    {
        private:
            std::vector<std::pair<Tm,Tc>> __data;
            mutable basic_polynomial_status __status;
            const compare* __comp=&init_comp;
        public:
            typedef compare compare_type;
            typedef Tm momomial_type;
            typedef Tc coeff_type;
            typedef Tm first_type;
            typedef Tc second_type;

            typedef typename std::vector<std::pair<Tm,Tc>>::iterator iterator;
            typedef typename std::vector<std::pair<Tm,Tc>>::const_iterator  const_iterator ;

            static compare init_comp;
            // static const Tc Tc_zero;
            basic_polynomial():__data(){}
            basic_polynomial(const compare *p):__data(),__comp(p){}
            basic_polynomial(const basic_polynomial &p)
            :__data(p.__data),__status(p.__status),__comp(p.__comp)
            {}
            basic_polynomial(basic_polynomial &&p)
            :__data(std::move(p.__data)),__status(std::move(p.__status)),__comp(p.__comp)
            {}
            basic_polynomial( std::initializer_list<std::pair<Tm,Tc>> init)
            :__data(init)
            {
                this->normalization();
            }
            basic_polynomial(const std::vector<std::pair<Tm,Tc>> & v)
            :__data(v)
            {
                this->normalization();
            }
            basic_polynomial(std::vector<std::pair<Tm,Tc>> && v)
            :__data(std::move(v))
            {
                this->normalization();
            }
            // basic_polynomial(const Tm& m)
            // :__data({{m,1}})
            // {}
            // basic_polynomial(Tm&& m)
            // :__data({std::pair<Tm,Tc>(std::move(m),1)})
            // {}
            // basic_polynomial(const Tc & c)
            // :__data({{Tm(),c}})
            // {}
            // basic_polynomial(Tc && c)
            // :__data({{Tm(),c}})
            // {}

            inline basic_polynomial & operator=(const basic_polynomial & p)
            {
                if (&p==this)
                    return *this;
                this->__data=p.__data;
                this->__status=p.__status;
                this->__comp=p.__comp;
                return *this;
            }

            inline basic_polynomial & operator=(basic_polynomial && p)
            {
                if (&p==this)
                    return *this;
                this->__data=std::move(p.__data);
                this->__status=std::move(p.__status);
                this->__comp=p.__comp;
                return *this;
            }

            inline basic_polynomial & operator=( std::initializer_list<std::pair<Tm,Tc>> init)
            {
                this->__data=init;
                this->__status.clear();
                this->normalization();
                return *this;
            }

            inline basic_polynomial & operator=(const std::vector<std::pair<Tm,Tc>> & init)
            {
                this->__data=init;
                this->__status.clear();
                this->normalization();
                return *this;
            }

            inline basic_polynomial & operator=(std::vector<std::pair<Tm,Tc>> && init)
            {
                this->__data=std::move(init);
                this->__status.clear();
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
            constexpr const compare & comp() const {return *(this->__comp);}
            //constexpr compare & comp() {return *(this->__comp);}
            inline bool comp(const variable & a,const variable &b)const {return (*this->__comp)(a,b);}
            constexpr void comp(const compare * c){this->__comp=c;}
            constexpr const compare * comp_ptr() const {return this->__comp;}

            constexpr void shrink_to_fit(){this->__data.shrink_to_fit();}
            constexpr std::size_t capacity(){return this->__data.capacity();}
            inline void clear() 
            {
                this->__data.clear();
                this->__status.clear();
            }

            inline void swap(basic_polynomial &p)
            {
                if (&p==this)
                    return void();
                this->__data.swap(p.__data);
                this->__status.swap(p.__swap);
                std::swap(this->__comp,p.__comp);
            }

            inline const std::list<std::pair<variable,int64_t>> & variables()  const
            {
                if (!this->__status.is_var)
                {
                    this->__status.is_var=true;
                    this->__status.variables=get_variables(*this);
                }
                return this->__status.variables;
            }

            inline int64_t degree() const
            {
                if (!this->__status.is_deg)
                {
                    this->__status.is_deg=true;
                    this->__status.deg=get_deg(*this);
                }
                return this->__status.deg;
            }

            inline bool is_normal()  const 
            {
                this->__status.clear();
                return pair_vec_first_normal_check(this->__data,this->comp()) && pair_vec_normal_check(this->__data,this->comp());
                // if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //     return false
                // for (auto && i:*this)
                //     if (!i.first.is_normal())
                //         return false;
                // return true;
            
            }

            inline void normalization()
            {
                this->__status.clear();
                // for (auto && i:*this)
                //     i.first.normalization());
                pair_vec_first_normalization(this->__data,this->comp());
                auto tmp_size=pair_vec_normalization(this->__data,pair_compare<Tm,Tc,compare>(this->__comp));
                this->resize(tmp_size);
                __auto_shrink(this->__data);

            }

            inline basic_polynomial operator+  (const basic_polynomial &p)const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
                // #ifdef DEBUG
                //     if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //         throw std::invalid_argument("Left basic_polynomial is not normal.");
                //     if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                //         throw std::invalid_argument("Right basic_polynomial is not normal.(In left comparation.)");
                // #endif
                basic_polynomial new_p;
                new_p.__comp=this->__comp;
                pair_vec_add(new_p.__data,this->__data,p.__data,this->comp());
                __auto_shrink(new_p.__data);
                return new_p;
            }

            // constexpr basic_polynomial operator+(const Tc & c) const
            // {
            //     return (*this)+basic_polynomial(c);
            // }

            // constexpr basic_polynomial operator-(const Tc & c) const
            // {
            //     return (*this)-basic_polynomial(c);
            // }


            inline basic_polynomial operator-  (const basic_polynomial &p)const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
                // #ifdef DEBUG
                //     if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //         throw std::invalid_argument("Left basic_polynomial is not normal.");
                //     if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                //         throw std::invalid_argument("Right basic_polynomial is not normal.(In left comparation.)");
                // #endif
                basic_polynomial new_p;
                new_p.__comp=this->__comp;
                pair_vec_sub(new_p.__data,this->__data,p.__data,this->comp());
                __auto_shrink(new_p.__data);
                return new_p;
            } 

            inline basic_polynomial operator+  () const
            {
                return *this;
            }   

            inline basic_polynomial operator- () const
            {
                basic_polynomial new_p;
                new_p.__comp=this->__comp;
                new_p.__data=pair_vec_negate(this->__data);
                return new_p;
            } 

            inline basic_polynomial operator*(const basic_polynomial & p) const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
                // #ifdef DEBUG
                //     if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //         throw std::invalid_argument("Left basic_polynomial is not normal.");
                //     if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                //         throw std::invalid_argument("Right basic_polynomial is not normal.(In left comparation.)");
                // #endif
                basic_polynomial new_p;
                new_p.__comp=this->__comp;
                pair_vec_multiplies(new_p.__data,this->__data,p.__data,this->comp());    
                __auto_shrink(new_p.__data);
                return new_p;
            }
            
            inline basic_polynomial operator/(const basic_polynomial & p) const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
                basic_polynomial new_p;
                new_p.__comp=this->__comp;
                pair_vec_div(new_p.__data,this->__data,p.__data,this->comp());    
                __auto_shrink(new_p.__data);
                return new_p;
            }

            inline basic_polynomial operator*(const Tc & p) const
            {
                basic_polynomial new_p;
                new_p.__comp=this->__comp;
                if (zore_check<Tc>()(p) || this->empty()) 
                    return new_p;
                new_p=*this;
                for (auto &i:new_p)
                    i.second*=p;
                return new_p;
            }

            inline basic_polynomial operator*(const Tm & m) const
            {
                basic_polynomial new_p;
                if (zore_check<Tm>()(m)) 
                {
                    new_p=*this;
                    return new_p;
                }
                new_p.__comp=this->__comp;
                for (auto &i:*this)
                {
                    new_p.push_back({i.first*m,i.second});
                }
                return new_p;
            }

            // inline basic_polynomial  power(unsigned i) const
            // {
            //     basic_polynomial newp;
            //     newp.__comp=this->__comp;
            //     if (this->empty() )
            //         return newp; 
            //     newp.push_back({Tm(this->__comp),1});
            //     if (i==0)
            //     {
            //         return newp;
            //     }
            //     basic_polynomial p(*this);
            //     while (i!=0)
            //     {
            //         //std::cout<<"power:"<<p<<std::endl;
            //         if (i%2!=0)
            //             newp=newp*p;
            //         i>>=1;
            //         if (i!=0)
            //             p=p*p;
            //     }
            //     return newp;
            // }


            inline bool operator> (const basic_polynomial &p) const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
               
                // #ifdef DEBUG
                //     if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //         throw std::invalid_argument("Left basic_polynomial is not normal.");
                //     if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                //         throw std::invalid_argument("Right basic_polynomial is not normal.(In left comparation.)");
                // #endif
                return pair_vec_comp(this->__data,p.__data,this->comp());
            }

            inline bool operator< (const basic_polynomial &p) const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
                // #ifdef DEBUG
                //     if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //         throw std::invalid_argument("Left basic_polynomial is not normal.");
                //     if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                //         throw std::invalid_argument("Right basic_polynomial is not normal.(In left comparation.)");
                // #endif
                return pair_vec_comp(p.__data,this->__data,this->comp());
            }

            inline bool operator== (const basic_polynomial &p) const
            {
                assert(this->__comp==p.__comp || this->comp()==p.comp());
                // #ifdef DEBUG
                //     if (!pair_vec_normal_check(this->begin(),this->end(),this->comp))
                //         throw std::invalid_argument("Left basic_polynomial is not normal.");
                //     if (!pair_vec_normal_check(p.begin(),p.end(),this->comp))
                //         throw std::invalid_argument("Right basic_polynomial is not normal.(In left comparation.)");
                // #endif
                return pair_vec_equal_to(p.__data,this->__data);
            }

            inline bool operator<=(const basic_polynomial &p) const
            {
                return !(*this>p);
            }

            inline bool operator>=(const basic_polynomial &p) const
            {
                return !(*this<p);
            }

            std::string str() const
            {
                if (this->size()==0)
                    return "0";
                std::ostringstream ss;
                ss<<(*this);
                return ss.str();
            }

            inline const_iterator find(const Tm & m) const
            {
                // std::function<bool(const std::pair<Tm,Tc> &,const std::pair<Tm,Tc> &)> tmp_comp=
                //     [&](std::pair<Tm,Tc> a,std::pair<Tm,Tc> b){return (this->comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<Tm,Tc>(m,0),pair_compare<Tm,Tc,compare>(this->__comp));
            }

            inline iterator find(const Tm & m)
            {
                // std::function<bool(const std::pair<Tm,Tc> &,const std::pair<Tm,Tc> &)> tmp_comp=
                //     [&](std::pair<Tm,Tc> a,std::pair<Tm,Tc> b){return (this->comp(a.first,b.first));};
                return pair_vec_find_first(this->begin(),this->end(),std::pair<Tm,Tc>(m,0),pair_compare<Tm,Tc,compare>(this->__comp));
            }

            // inline const Tc & coef(const Tm & m) const
            // {
            //     auto l=this->find(m);
            //     if (l!=this->end())
            //         return l->second;
            //     else
            //     {
            //         //Tc tmp_Tc=0;
            //         return Tc_zero;//tmp_Tc;
            //     }
                    
            // }

            friend std::ostream& operator<<  (std::ostream& stream, const basic_polynomial& v) 
            {
                if (v.size()==0)
                    return stream<<'0';
                bool is_print=false;
                for(auto i=v.begin();i!=v.end();++i)
                {
                    if (!zore_check<Tm>()(i->first) //&& !zore_check<Tc>()(i->second)
                        ){
                        if (is_print && i->second>=0) //greater
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
                    if (zore_check<Tm>()(i->first) //&& !zore_check<Tc>()(i->second)
                        )
                        if (is_print && i->second>=0) //greater
                            stream<<"+";
                        is_print=true;
                        stream<<i->second;
                }
                return stream;
            }
              
    };

    template<class Tm,class Tc,class compare>  compare  basic_polynomial<Tm, Tc, compare>::init_comp=compare();
    // template<class Tm,class Tc,class compare> const Tc  basic_polynomial<Tm, Tc, compare>::Tc_zero=0;
    
    template<class Tm,class Tc,class compare>
    struct zore_check<basic_polynomial<Tm, Tc, compare>>: public std::unary_function<basic_polynomial<Tm, Tc, compare>, bool>
    {
        constexpr bool operator()(const basic_polynomial<Tm, Tc, compare> & m)
        {
            return m.empty();
        } 
    };
    template<class T1,class T2,class T3,class T4,class comp>
    inline void polynomial_div(basic_polynomial<T1,T2,comp> & p,const basic_polynomial<T1,T3,comp> & p1,const basic_polynomial<T1,T4,comp> & p2)
    {
        assert(p1.comp_ptr()==p2.comp_ptr() || p1.comp()==p2.comp());
        p.comp(p1.comp_ptr());
        pair_vec_div(p.data(),p1.data(),p2.data(),p1.comp());    
        __auto_shrink(p.data());
    }
    
    // template<class Tc,class Tm,class compare>
    // void swap(basic_polynomial<Tm,Tc,compare> & v1,basic_polynomial<Tm,Tc,compare> & v2)
    // {
    //     v1.swap(v2);
    // }


}
namespace std{
    template<class Tm,class Tc,class compare>
    struct hash<clpoly::basic_polynomial<Tm,Tc,compare>>
    {
        std::size_t operator()(const clpoly::basic_polynomial<Tm,Tc,compare> & m) const
        {
            return (hash<std::string>()(m.str()));
        }
    };
}
#endif