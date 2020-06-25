/*
Module Name:
    monomial.hh
Abstract:
    定义类：monomial
    monomial有关得操作实现
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
#include <cassert>
#include "monomial_order.hh"
#include "basic_monomial.hh"
namespace clpoly
{
    using monomial=basic_monomial<grlex>;
    
    template<class compare>
    struct zore_check<basic_monomial<compare>>: public std::unary_function<basic_monomial<compare>, bool>
    {
        constexpr bool operator()(const basic_monomial<compare> & m)
        {
            return m.empty();
        } 
    };
    
    template <class T,class compare>
    inline void pair_vec_first_normalization(std::vector<std::pair<basic_monomial<compare>,T>> & v,const compare * comp)
    {   
        for(auto &i:v)
        {
            // if (i.first.comp_ptr()!=comp)
            // {
                i.first.comp(comp);
                i.first.normalization();
            // }
        }
    }
    
    template <class T,class compare>
    inline bool pair_vec_first_normal_check(std::vector<std::pair<basic_monomial<compare>,T>> & v,const  compare * comp)
    {
        for(auto &i:v)
            if (i.first.comp_ptr()!=comp || i.is_normal())
            {
                return false;
            }
        return true;
    }
    
    template<class compare>
    void __mono_mult__(basic_monomial<compare> & op,const basic_monomial<compare>  & op1,const basic_monomial<compare>  & op2)
    {
        // if(&(op.data())==&(op1.data()) || &(op.data())==&(op2.data()))
        //     op=op1*op2;
        //std::cout<<"__mono_mult___1";
        op.comp(op1.comp_ptr());
        pair_vec_add(op.data(),op1.data(),op2.data(),op1.comp());
        op.deg()=op1.deg()+op2.deg();
        //std::cout<<op<<"="<<op1<<"*"<<op2<<std::endl;
    }
    
    template<class compare>
    bool is_divexact(basic_monomial<compare> & op,const basic_monomial<compare>  & op1,const basic_monomial<compare>  & op2)
    {
        op.comp(op1.comp_ptr());
        pair_vec_sub(op.data(),op1.data(),op2.data(),op1.comp());
        op.deg()=op1.deg()-op2.deg();
        if (op.deg()<0) return false;
        for(auto &i : op)
            if (i.second < 0)
                return false;
        return true;
    }
    
    template<class Tc,class compare>
    void __pair_vec_variables(const std::vector<std::pair<basic_monomial<compare>,Tc>> & p,std::list<std::pair<variable,int64_t>>& l)
    {
        typename std::list<std::pair<variable,int64_t>>::iterator l_ptr;
        typename basic_monomial<compare>::const_iterator m_ptr;

        for (const auto & i:p)
        {
            if (!i.first.empty())
            {
                m_ptr=i.first.begin();
                if (!l.empty())
                {
                    l_ptr=l.begin();
                    while(m_ptr!=i.first.end())
                    {
                        while (l_ptr!=l.end() && i.first.comp(l_ptr->first,m_ptr->first))
                        {
                            ++l_ptr;
                        }
                        if (l_ptr==l.end()) break;
                        if (i.first.comp(m_ptr->first,l_ptr->first))
                        {
                            l.insert(l_ptr,*m_ptr);
                        }
                        else 
                        {
                            if (m_ptr->second>l_ptr->second)
                                l_ptr->second=m_ptr->second;
                            ++l_ptr;
                        }    
                        ++m_ptr;
                    }
                }
                for(;m_ptr!=i.first.end();l.push_back(*(m_ptr++)));  
            }
        }
    }
    template<class Tc,class compare>
    void __pair_vec_variables(const std::vector<std::pair<basic_monomial<compare>,Tc>> & p,std::list<variable>& l)
    {
        typename std::list<variable>::iterator l_ptr;
        typename basic_monomial<compare>::const_iterator m_ptr;

        for (const auto & i:p)
        {
            if (!i.first.empty())
            {
                m_ptr=i.first.begin();
                if (!l.empty())
                {
                    l_ptr=l.begin();
                    while(m_ptr!=i.first.end())
                    {
                        while (l_ptr!=l.end() && i.first.comp(*l_ptr,m_ptr->first))
                        {
                            ++l_ptr;
                        }
                        if (l_ptr==l.end()) break;
                        if (i.first.comp(m_ptr->first,*l_ptr))
                        {
                            l.insert(l_ptr,m_ptr->first);
                        }
                        else 
                        {
                            ++l_ptr;
                        }    
                        ++m_ptr;
                    }
                }
                for(;m_ptr!=i.first.end();l.push_back((m_ptr++)->first));  
            }
        }
    }
    
  
    bool __pair_vec_multiplies_compression_b=false;
    template <class T2,class T3,class T4,class compare>
    bool __pair_vec_multiplies_compression
    (
        std::vector<std::pair<basic_monomial<compare>,T2>>& new_v,
        const std::vector<std::pair<basic_monomial<compare>,T3>> & v1_,
        const std::vector<std::pair<basic_monomial<compare>,T4>> & v2_,
        const compare & comp
    )
    {
        if (v1_.empty() || v2_.empty()) return false;
        if (!__pair_vec_multiplies_compression_b) return false; 
        std::list<variable> vars;
        std::cout<<"压缩乘法开始\n";
        if (!__is_monomial_multiplies_compression(v1_,v2_,comp,vars)) return false;
        std::cout<<"准备完成\n";
        //std::cout<<"vars:"<<vars.size()<<std::endl;
        std::greater<uint64_t> gcomp;
        // std::vector<std::pair<uint64_t,const T3*>> v1_c;
        std::vector<std::pair<uint64_t,const T4*>> v2_c;
        v2_c.reserve(v2_.size());
        // std::vector<std::pair<uint64_t,T2*>> new_v_c;
        new_v.reserve(v1_.size()*v2_.size());
        const compare * comp_ptr=v1_.begin()->first.comp_ptr();
        VHC<uint64_t,std::pair<uint64_t,const T3*>,typename std::vector<std::pair<uint64_t,const T4*>>::const_iterator > **heap=
            new VHC<uint64_t,std::pair<uint64_t,const T3*>,typename std::vector<std::pair<uint64_t,const T4*>>::const_iterator >*[v1_.size()];
        VHC<uint64_t,std::pair<uint64_t,const T3*>,typename std::vector<std::pair<uint64_t,const T4*>>::const_iterator > *node=
            new VHC<uint64_t,std::pair<uint64_t,const T3*>,typename std::vector<std::pair<uint64_t,const T4*>>::const_iterator >[v1_.size()];    
        VHC<uint64_t,std::pair<uint64_t,const T3*>,typename std::vector<std::pair<uint64_t,const T4*>>::const_iterator > **lin=
            new VHC<uint64_t,std::pair<uint64_t,const T3*>,typename std::vector<std::pair<uint64_t,const T4*>>::const_iterator >*[v1_.size()];
        std::size_t reset=1;
        std::size_t lin_size=0;
        auto v1_ptr=v1_.begin();
        for (auto &i:v2_)
        {
            v2_c.push_back({_monomial_compression(i.first,vars),&(i.second)});
        }
        auto v2_begin=v2_c.begin();
        auto v2_end=v2_c.end();
        for(std::size_t i=0;i!=v1_.size();++i,++v1_ptr)
        {
            heap[i]=node+i;
            node[i].v2_ptr=v2_begin;
            node[i].v1_ptr.first=_monomial_compression(v1_ptr->first,vars);
            node[i].v1_ptr.second=&v1_ptr->second;
            node[i].mono=v2_begin->first+node[i].v1_ptr.first;
            //node[i].next=nullptr;
        } 
        heap[0]->next=nullptr;
        std::size_t heap_size=1;
        std::size_t i, j, s,i1;
        uint64_t m;
        basic_monomial<compare> m1;
        T2 k;
        while (heap_size!=0)
        {
            set_zero(k);
            m=heap[0]->mono;
            // _monomial_decompression(m,m1,vars,comp_ptr);
            // std::cout<<m1<<std::endl;
            do{
                while(heap[0]!=nullptr){
                    addmul(k,*(heap[0]->v1_ptr.second),*(heap[0]->v2_ptr->second));
                    if (heap[0]->v2_ptr==v2_begin&& reset!=v1_.size()){
                        lin[lin_size++]=heap[reset++];
                    }
                    if (++(heap[0]->v2_ptr)!=v2_end){
                        lin[lin_size++]=heap[0];
                        heap[0]->mono=heap[0]->v2_ptr->first+heap[0]->v1_ptr.first;
                        heap[0]=heap[0]->next;  
                    }
                    else{
                        heap[0]=heap[0]->next;
                    }
                }
                VHC_extract(heap,heap_size,gcomp);
                //end:
            }while (heap_size>0 && heap[0]->mono==m);  //equal_to
            if (!zore_check<T2>()(k)){
                _monomial_decompression(m,m1,vars,comp_ptr);
                new_v.push_back({std::move(m1),std::move(k)});
            }
            while(lin_size>0)
                VHC_insert(heap,heap_size,lin[--lin_size],gcomp);
        }
        delete [] heap;
        delete [] lin;
        delete [] node;
        std::cout<<"压缩乘法结束\n";
        return true;
    }
}

#endif