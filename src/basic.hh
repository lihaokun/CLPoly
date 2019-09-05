/*
Module Name:
    basic.hh
Abstract:
    这里定义一些基础元素
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_BASIC_HH
#define CLPOLY_BASIC_HH
#include <vector>
#include <algorithm>
#include <functional>
#define shrink_bound 0.7
#define shrink_bound_abs 100

namespace clpoly{
    template<class T>
    inline bool zore_check(const T & c)
    {
        return c==0;
    }
    template <class T1,class T2>
    inline T1 plus(const T1 & op1,const T2&op2)
    {
        return op1+op2;
    }
    template <class T1,class T2>
    inline T1 multiplies(const T1 & op1,const T2&op2)
    {
        return op1*op2;
    }
    template<class T1,class T2>
    inline void add_assignment(T1 & x,const T2 & y)
    {
        x+=y;
    } 
    template<class T1,class T2>
    inline void sub_assignment(T1 & x,const T2 & y)
    {
        x-=y;
    } 
    template <class T1,class T2,class T3>
    inline void addmul(T1 & op,const T2&op1,const T3&op2)
    {
        op+=op1*op2;
    }
    template <class T1,class T2,class T3>
    inline void submul(T1 & op,const T2&op1,const T3&op2)
    {
        op-=op1*op2;
    }
    template<class T>
    inline T negate (const T & m)
    {
        return -m;
    }
    template<class T1,class T2>
    inline bool greater(const T1 & m1,const T2 & m2)
    {
        return m1>m2;
    }
    template<class T1,class T2>
    inline bool less(const T1 & m1,const T2 & m2)
    {
        return m1<m2;
    }
    template<class T1,class T2>
    inline bool  equal_to(const T1 & m1,const T2 & m2)
    {
        return m1==m2;
    }
    template<class T>
    inline void set_zero(T & x)
    {
        x=0;
    } 
    template <class T,class compare>
    inline bool pair_vec_normal_check(const T & b,const T & e,const compare & comp )
    {
        T t2=b;++t2;
        T t1=b;
        if (zore_check(t1->second))
            return false;
        while (t2!=e)
        {
            if (!comp(t1->first,t2->first) || zore_check(t2->second))
                return false;
            t2++;
        }
        return true;
    }
    template <class T,class compare>
    std::size_t pair_vec_normalization(const T & b,const T & e,const compare & comp )
    {
        if (!std::is_sorted(b,e,comp))
            std::sort(b,e,comp);
        T t2=b;
        while (t2!=e && zore_check(t2->second))
            ++t2;
        if (t2==e)
            return 0;
        T t1=b;
        //std::cout<<t2->first<<":"<<t2->second<<std::endl;
        if (t1!=t2)
            *t1=std::move(*t2);
        ++t2;
        std::size_t tmp_size=1;
        for(;t2!=e;++t2)
            //std::cout<<t2->first<<":"<<t2->second<<std::endl;
            if (!zore_check(t2->second))
                if ( equal_to(t1->first,t2->first))
                    add_assignment(t1->second,t2->second);
                else
                {
                    ++t1;
                    ++tmp_size;
                    if (t1!=t2)
                        *t1=std::move(*t2);
                }
        return tmp_size;
    }
    template <class T>
    T pair_vec_negate(const T & v1)
    {
        T new_T;
        new_T.reserve(v1.size());
        for (auto &&i:v1)
        {
            new_T.push_back({i.first,negate(i.second)});
        }
        return new_T;
    }
    template <class T,class compare>
    T pair_vec_add(const T & v1,const T & v2,const compare & comp)
    {
        if (v1.empty())
            return v2;
        if (v2.empty())
            return v1;
        T v;
        v.reserve(std::max(v1.size(),v2.size()));
        auto v1_ptr=v1.begin();
        auto v2_ptr=v2.begin();
        auto v1_end=v1.end();
        auto v2_end=v2.end();
        
        while (v1_ptr!=v1_end && v2_ptr!=v2_end)
        {
            if (comp(v2_ptr->first,v1_ptr->first))
                v.push_back(*(v2_ptr++));
            else
            {
                v.push_back(*(v1_ptr++));
                if ( equal_to(v2_ptr->first,v.back().first))
                {
                    add_assignment(v.back().second,(v2_ptr++)->second);
                    if (zore_check(v.back().second))
                        v.pop_back();

                }
            }
        }
    
        if (v2_ptr!=v2_end)
        {
            v1_ptr=v2_ptr;
            v1_end=v2_end;
        }
        v.reserve(v.size()+(v1_end-v1_ptr)+1);
        while (v1_ptr!=v1_end)
        {
            v.push_back(*(v1_ptr++));
        }
        return v;
    }
    template <class T,class compare>
    T pair_vec_sub(const T & v1,const T & v2,const compare & comp)
    {
        if (v1.empty())
            return pair_vec_negate(v2);
        if (v2.empty())
            return v1;
        T v;
        v.reserve(std::max(v1.size(),v2.size()));
        auto v1_ptr=v1.begin();
        auto v2_ptr=v2.begin();
        auto v1_end=v1.end();
        auto v2_end=v2.end();
        
        while (v1_ptr!=v1_end && v2_ptr!=v2_end)
        {
            if (comp(v2_ptr->first,v1_ptr->first))
            {
                v.push_back({v2_ptr->first,negate(v2_ptr->second)});
                ++v2_ptr;
            }
            else{
                v.push_back(*(v1_ptr++));
                if ( equal_to(v2_ptr->first,v.back().first))
                {
                    sub_assignment(v.back().second,(v2_ptr++)->second);
                    if (zore_check(v.back().second))
                        v.pop_back();
                }
            }
        }
    
        if (v2_ptr!=v2_end)
        {
            v.reserve(v.size()+(v2_end-v2_ptr)+1);
            
            while (v2_ptr!=v2_end)
            {
                v.push_back({v2_ptr->first,negate(v2_ptr->second)});
                ++v2_ptr;
            }
        }
        else
        {
            v.reserve(v.size()+(v1_end-v1_ptr)+1);
            while (v1_ptr!=v1_end)
            {
                v.push_back(*(v1_ptr++));
            }
        }
        return v;
    }
    template <class T,class compare1,class compare2>
    bool pair_vec_comp(const T & v1,const T & v2,const compare1 & comp1,const compare2 & comp2)
    {
        auto v1_ptr=v1.begin();
        auto v1_end=v1.end();
        auto v2_ptr=v2.begin();
        auto v2_end=v2.end();
        for(;(v1_ptr!=v1_end && v2_ptr!=v2_end);++v1_ptr,++v2_ptr)
        {
            if (!equal_to(v1_ptr->first,v2_ptr->first))
                return comp1(v1_ptr->first,v2_ptr->first);
            if (!equal_to(v1_ptr->second,v2_ptr->second))
                return comp2(v1_ptr->second,v2_ptr->second);
        }
        if (v1_ptr==v1_end && v2_ptr!=v2_end)
            return true;
        else
            return false;
    }
    template <class T,class compare1,class compare2>
    bool pair_vec_equal_to(const T & v1,const T & v2,const compare1 & comp1,const compare2 & comp2)
    {
        if (v1.size()!=v2.size())
            return false;
        auto v1_ptr=v1.begin();
        auto v1_end=v1.end();
        auto v2_ptr=v2.begin();
        auto v2_end=v2.end();
        for(;(v1_ptr!=v1_end && v2_ptr!=v2_end);++v1_ptr,++v2_ptr)
        {
            if (!equal_to(v1_ptr->first,v2_ptr->first) || !equal_to(v1_ptr->second,v2_ptr->second))
                return false;
        }
        return true;
    }
    template<class T,class ptr1,class ptr2>
    struct VHC{
        T mono;
        ptr1 v1_ptr;
        ptr2 v2_ptr;
        VHC<T,ptr1,ptr2>* next;
    };
    template<class T,class ptr1,class ptr2>
    inline VHC<T,ptr1,ptr2>* VHC_init(VHC<T,ptr1,ptr2>* H,ptr1 v1_ptr,ptr2 v2_ptr)
    {
        H->mono=multiplies(v1_ptr->first,v2_ptr->first);
        H->v1_ptr=v1_ptr;
        H->v2_ptr=v2_ptr;
        H->next=nullptr;
        return H;
    }
    template<class T,class ptr1,class ptr2,class compare>
    inline void VHC_extract(VHC<T,ptr1,ptr2>** heap,std::size_t & heap_size,const compare & comp)
    {
        std::size_t i,j,s;
        s = --heap_size;
        for (i=0, j=1; j < s; i=j, j=(j<<1)+1) {
            j = (comp(heap[j]->mono , heap[j+1]->mono)) ? j : j+1;
            if (comp(heap[j]->mono,heap[s]->mono))
                heap[i] = heap[j];
            else{
                break;  
            }
        }
        heap[i] =heap[s];
        
    }
    template<class T,class ptr1,class ptr2,class compare>
    inline void VHC_insert(VHC<T,ptr1,ptr2>** heap,std::size_t & heap_size,VHC<T,ptr1,ptr2>* newVHC,const compare & comp)
    {
        std::size_t i,j,i1;
        if (heap_size==0)
        {
            heap_size++;
            heap[0]=newVHC;
            newVHC->next=nullptr;
            //return void();
            }
            else
            if (equal_to(newVHC->mono,heap[0]->mono)){
                newVHC->next=heap[0];
                heap[0]=newVHC;
                //return void();
            }
            else
                if (comp(newVHC->mono,heap[0]->mono)){
                for (i=heap_size++, j=(i-1)>>1; i>0; heap[i]=heap[j], i=j, j=(j-1)>>1);
                heap[0]=newVHC;newVHC->next=nullptr;
                //return void();
                }
                else{
                    for(i1=(heap_size-1)>>1;comp(newVHC->mono ,heap[i1]->mono);i1=(i1-1)>>1);
                    if (newVHC->mono == heap[i1]->mono)
                    {
                        newVHC->next=heap[i1];
                        heap[i1]=newVHC;
                    }
                    else{
                        for (i=heap_size++, j=(i-1)>>1; j!=i1; heap[i]=heap[j], i=j, j=(j-1)>>1);
                        heap[i]=newVHC;newVHC->next=nullptr;
                    }
                }
    }

    template <class T1,class T2,class compare>
    std::vector<std::pair<T1,T2>> pair_vec_multiplies
    (
        const std::vector<std::pair<T1,T2>> & v1_,
        const std::vector<std::pair<T1,T2>> & v2_,
        const compare & comp
    )
    {
        std::vector<std::pair<T1,T2>> new_v;
        const std::vector<std::pair<T1,T2>> * v1，*v2;
        if (v1_.size()>v2_.size())
        {
            v1=&v2_;v2=&v1_;
        }
        else
        {
            v1=&v1_;v2=&v2_;
        }
        if (v1->size()==0)
            return new_v;
        if (v1->size()==1)
        {
            new_v.reserve(v2->size());
            for(auto &&i:*v2)
            {
                new_v.push_back(
                    multiplies((*v1)[0].first,i.first),
                    multiplies((*v1)[0].second,i.second));
                
            }
            return new_v;
        }
        VHC<T1,std::pair<T1,T2>*,std::pair<T1,Type>*>* heap[v1->size()];
        VHC<T1,std::pair<T1,T2>*,std::pair<T1,Type>*>* lin[v1->size()];
        std::size_t reset=1;
        VHC<T1,std::pair<T1,Type>*,std::pair<T1,Type>*> *lout;
        std::size_t lin_size=0;
        auto v1_ptr=v1->begin();
        auto v2_begin=v2->begin();
        auto v2_end=v2->end();
        for(std::size_t i=0;i!=v1->size();++i)
            heap[i]=VHC_init(new VHC<T1,std::pair<T1,T2>*,std::pair<T1,T2>*>,v1_ptr++,v2_begin);
        std::size_t heap_size=1;
        std::size_t i, j, s,i1;
        T1 m;
        T2 k;
        new_v.reserve(v1->size()+v2->size());
        while (heap_size!=0)
        {
            set_zero(k);
            m=std::move(heap[0]->mono);
            do{
                while(heap[0]!=nullptr){
                    addmul(k,heap[0]->v1_ptr->second,heap[0]->v2_ptr->second);
                    if (heap[0]->v2_ptr==v2_begin&& reset!=p1.size()){
                        lin[lin_size++]=heap[reset++];
                    }
                    if (++(heap[0]->v2_ptr)!=v2_end){
                        lin[lin_size++]=heap[0];
                        heap[0]->mono=multiplies(heap[0]->v1_ptr->first,heap[0]->v2_ptr->first);
                        heap[0]=heap[0]->next;  
                    }
                    else{
                        lout=heap[0];
                        delete lout;
                        heap[0]=heap[0]->next;
                    }
                }
                VHC_extract(heap,heap_size);
                //end:
            }while (heap_size>0 && equal_to(heap[0]->mono,m));
            if (k!=0){
                new_v.push_back({std::move(m),std::move(k)});
            }
            while(lin_size>0)
                VHC_insert(heap,heap_size,lin[--lin_size]);
        }
        if (new_v.size()*1.0/new_v.capacity()<shrink_bound && new_v.size()>shrink_bound_abs)
        {
            std::cont<<"shrink_to_fit";
            new_v.shrink_to_fit();
        }
        return new_v;  
    }
}
#endif