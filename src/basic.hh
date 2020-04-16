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
#include <iostream>
#define CLPOLY_shrink_bound 0.5
#define CLPOLY_shrink_bound_abs 100

namespace clpoly{
    template<class T>
    struct zore_check: public std::unary_function<T, bool>
    {
        constexpr bool operator()(const T & c)
        {
            return c==0;
        } 
    };

    template <class T1,class T2,class T3>
    constexpr void addmul(T1 & op,const T2&op1,const T3&op2)
    {
        op+=op1*op2;
    }
    template <class T1,class T2,class T3>
    constexpr void submul(T1 & op,const T2&op1,const T3&op2) 
    {
        op-=op1*op2;
    }
    using std::less;
    using std::greater;
    template <class T1>
    constexpr void set_zero(T1& op)
    {
        op=0;
    }
    template <class T1,class T2,class compare>
    inline bool pair_vec_normal_check(const std::vector<std::pair<T1,T2>> & v,const compare & comp )
    {
        auto t2=v.begin();++t2;
        auto t1=v.begin();
        const auto & e=v.end();
        if (zore_check<T2>()(t1->second))
            return false;
        while (t2!=e)
        {
            if (!comp(t1->first,t2->first) || zore_check<T2>()(t2->second))
                return false;
            ++t2;++t1;
        }
        return true;
    }
    template <class T1,class T2,class compare>
    std::size_t pair_vec_normalization(std::vector<std::pair<T1,T2>> & v,const compare & comp )
    {
        const auto &b=v.begin();
        const auto &e=v.end();
        if (!std::is_sorted(b,e,comp))
            std::sort(b,e,comp);
        auto t2=b;
        while (t2!=e && zore_check<T2>()(t2->second))
            ++t2;
        if (t2==e)
            return 0;
        auto t1=b;
        //std::cout<<t2->first<<":"<<t2->second<<std::endl;
        if (t1!=t2)
            *t1=std::move(*t2);
        ++t2;
        std::size_t tmp_size=1;
        for(;t2!=e;++t2)
            //std::cout<<t2->first<<":"<<t2->second<<std::endl;
            if (!zore_check<T2>()(t2->second))
                if (t1->first==t2->first) //equal_to
                    t1->second+=t2->second;//add_assignment
                else
                {
                    ++t1;
                    ++tmp_size;
                    if (t1!=t2)
                        *t1=std::move(*t2);
                }
        return tmp_size;
    }
    template <class T,class Tv,class compare>
    inline T pair_vec_find_first(const T & b,const T & e,const Tv & value,const compare & comp )
    {
        auto tmp_l = std::lower_bound(b, e, value, comp);
        return (tmp_l!=e && (value.first==tmp_l->first)) ? tmp_l : e;// equal_to
    }
    
    template <class T1,class T2>
    std::vector<std::pair<T1,T2>> pair_vec_negate(const std::vector<std::pair<T1,T2>> & v1)
    {
        std::vector<std::pair<T1,T2>> new_T;
        new_T.reserve(v1.size());
        for (auto &&i:v1)
        {
            new_T.push_back(std::pair<T1,T2>(i.first,-i.second));//negate
        }
        return new_T;
    }
    template <class T1,class T2,class compare>
    std::vector<std::pair<T1,T2>> pair_vec_add(
        const std::vector<std::pair<T1,T2>> & v1,const std::vector<std::pair<T1,T2>> & v2,const compare & comp)
    {
        if (v1.empty())
            return v2;
        if (v2.empty())
            return v1;
        std::vector<std::pair<T1,T2>> v;
        v.reserve(v1.size()+v2.size());
        auto v1_ptr=v1.begin();
        auto v2_ptr=v2.begin();
        const auto &v1_end=v1.end();
        const auto &v2_end=v2.end();
        T2 tmp;
        while (v1_ptr!=v1_end && v2_ptr!=v2_end)
        {
            if (comp(v2_ptr->first,v1_ptr->first))
                v.push_back(*(v2_ptr++));
            else
            {
                
                // v.push_back(*(v1_ptr++));
                // //std::cout<<v.back().first<<" ";
                // if ( equal_to<T1>()(v2_ptr->first,v.back().first))
                // {
                //     //std::cout<<"1 ";    
                //     add_assignment(v.back().second,(v2_ptr++)->second);
                //     if (zore_check<T2>()(v.back().second))
                //         v.pop_back();

                // }
                if ( v2_ptr->first==v1_ptr->first) //equal_to
                {
                    if (!zore_check<T2>()(tmp=v2_ptr->second+v1_ptr->second)) //plus
                        v.push_back(std::pair<T1,T2>(v1_ptr->first,tmp));
                    ++v2_ptr;++v1_ptr;
                }
                else
                    v.push_back(*(v1_ptr++)); 
            }
        }
    
        // if (v2_ptr!=v2_end)
        // {
        //     v1_ptr=v2_ptr;
        //     v1_end=v2_end;
        // }
        // //v.reserve(v.size()+(v1_end-v1_ptr));
        // while (v1_ptr!=v1_end)
        // {
        //     v.push_back(*(v1_ptr++));
        // }
        while (v1_ptr!=v1_end)
        {
            v.push_back(*(v1_ptr++));
        }
        while (v2_ptr!=v2_end)
        {
            v.push_back(*(v2_ptr++));
        }

        if (v.size()*1.0/v.capacity()<CLPOLY_shrink_bound && v.size()>CLPOLY_shrink_bound_abs)
        {
            //std::cout<<"shrink_to_fit";
            v.shrink_to_fit();
        }
        
        return v;
    }

    template <class T1,class T2,class compare>
    std::vector<std::pair<T1,T2>> pair_vec_sub(
        const std::vector<std::pair<T1,T2>> & v1,const std::vector<std::pair<T1,T2>> & v2,const compare & comp)
    {
        if (v1.empty())
            return pair_vec_negate(v2);
        if (v2.empty())
            return v1;
        std::vector<std::pair<T1,T2>> v;
        v.reserve(v1.size()+v2.size());
        auto v1_ptr=v1.begin();
        auto v2_ptr=v2.begin();
        const auto & v1_end=v1.end();
        const auto & v2_end=v2.end();
        T2 tmp;

        while (v1_ptr!=v1_end && v2_ptr!=v2_end)
        {
            if (comp(v2_ptr->first,v1_ptr->first))
            {
                v.push_back(std::pair<T1,T2>(v2_ptr->first,-v2_ptr->second));//negate
                ++v2_ptr;
            }
            else{
                // v.push_back(*(v1_ptr++));
                // if ( equal_to(v2_ptr->first,v.back().first))
                // {
                //     sub_assignment(v.back().second,(v2_ptr++)->second);
                //     if (zore_check(v.back().second))
                //         v.pop_back();
                // }
                if ( v2_ptr->first==v1_ptr->first) //equal_to
                {
                    if (!zore_check<T2>()(tmp=v1_ptr->second-v2_ptr->second)) //minus
                        v.push_back(std::pair<T1,T2>(v1_ptr->first,tmp));
                    ++v2_ptr;++v1_ptr;
                }
                else
                    v.push_back(*(v1_ptr++)); 

            }
        }
    
        if (v2_ptr!=v2_end)
        {
            //v.reserve(v.size()+(v2_end-v2_ptr)+1);
            
            while (v2_ptr!=v2_end)
            {
                v.push_back(std::pair<T1,T2>(v2_ptr->first,-v2_ptr->second));//negate
                ++v2_ptr;
            }
        }
        else
        {
            //v.reserve(v.size()+(v1_end-v1_ptr)+1);
            while (v1_ptr!=v1_end)
            {
                v.push_back(*(v1_ptr++));
            }
        }

        if (v.size()*1.0/v.capacity()<CLPOLY_shrink_bound && v.size()>CLPOLY_shrink_bound_abs)
        {
            //std::cout<<"shrink_to_fit";
            v.shrink_to_fit();
        }

        return v;
    }
    template <class T1,class T2,class compare>
    bool pair_vec_comp(const std::vector<std::pair<T1,T2>> & v1,const std::vector<std::pair<T1,T2>> & v2,const compare & comp)
    {
        auto v1_ptr=v1.begin();
        const auto & v1_end=v1.end();
        auto v2_ptr=v2.begin();
        const auto & v2_end=v2.end();
        for(;(v1_ptr!=v1_end && v2_ptr!=v2_end);++v1_ptr,++v2_ptr)
        {
            if (v1_ptr->first!=v2_ptr->first) //equal_to
                return comp(v1_ptr->first,v2_ptr->first);
            if (v1_ptr->second!=v2_ptr->second)//equal_to
                return v1_ptr->second>v2_ptr->second; //greater
        }
        if (v1_ptr==v1_end)
            return false;
        else
            return true;
    }
    template <class T>
    bool pair_vec_equal_to(const T & v1,const T & v2)
    {
        if (v1.size()!=v2.size())
            return false;
        auto v1_ptr=v1.begin();
        auto v1_end=v1.end();
        auto v2_ptr=v2.begin();
        auto v2_end=v2.end();
        for(;(v1_ptr!=v1_end && v2_ptr!=v2_end);++v1_ptr,++v2_ptr)
        {
            if (v1_ptr->first!=v2_ptr->first || v1_ptr->second!=v2_ptr->second) //equal_to
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
        H->mono=v1_ptr->first*v2_ptr->first; //multiplies
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
            if (newVHC->mono==heap[0]->mono){  //equal_to
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
        const std::vector<std::pair<T1,T2>> * v1,*v2;
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
                new_v.push_back({
                    (*v1)[0].first*i.first,     //multiplies
                    (*v1)[0].second*i.second}); //multiplies
                
            }
            return new_v;
        }
        VHC<T1,typename std::vector<std::pair<T1,T2>>::const_iterator,typename std::vector<std::pair<T1,T2>>::const_iterator > **heap=
            new VHC<T1,typename std::vector<std::pair<T1,T2>>::const_iterator,typename std::vector<std::pair<T1,T2>>::const_iterator >*[v1->size()];
        VHC<T1,typename std::vector<std::pair<T1,T2>>::const_iterator ,typename std::vector<std::pair<T1,T2>>::const_iterator > **lin=
            new VHC<T1,typename std::vector<std::pair<T1,T2>>::const_iterator,typename std::vector<std::pair<T1,T2>>::const_iterator >*[v1->size()];
        std::size_t reset=1;
        VHC<T1,typename std::vector<std::pair<T1,T2>>::const_iterator ,typename std::vector<std::pair<T1,T2>>::const_iterator > *lout;
        std::size_t lin_size=0;
        auto v1_ptr=v1->begin();
        auto v2_begin=v2->begin();
        auto v2_end=v2->end();
        for(std::size_t i=0;i!=v1->size();++i)
            heap[i]=VHC_init(new VHC<T1,typename std::vector<std::pair<T1,T2>>::const_iterator,typename std::vector<std::pair<T1,T2>>::const_iterator>,v1_ptr++,v2_begin);
        std::size_t heap_size=1;
        std::size_t i, j, s,i1;
        T1 m;
        T2 k;
        new_v.reserve(v1->size()*v2->size());
        while (heap_size!=0)
        {
            set_zero(k);
            m=std::move(heap[0]->mono);
            do{
                while(heap[0]!=nullptr){
                    addmul(k,heap[0]->v1_ptr->second,heap[0]->v2_ptr->second);
                    if (heap[0]->v2_ptr==v2_begin&& reset!=v1->size()){
                        lin[lin_size++]=heap[reset++];
                    }
                    if (++(heap[0]->v2_ptr)!=v2_end){
                        lin[lin_size++]=heap[0];
                        heap[0]->mono=heap[0]->v1_ptr->first*heap[0]->v2_ptr->first; //multiplies
                        heap[0]=heap[0]->next;  
                    }
                    else{
                        lout=heap[0];
                        delete lout;
                        heap[0]=heap[0]->next;
                    }
                }
                VHC_extract(heap,heap_size,comp);
                //end:
            }while (heap_size>0 && heap[0]->mono==m);  //equal_to
            if (k!=0){
                new_v.push_back({std::move(m),std::move(k)});
            }
            while(lin_size>0)
                VHC_insert(heap,heap_size,lin[--lin_size],comp);
        }
        if (new_v.size()*1.0/new_v.capacity()<CLPOLY_shrink_bound && new_v.size()>CLPOLY_shrink_bound_abs)
        {
            
            //std::cout<<"shrink_to_fit size:"<<new_v.size()<<" capacity:"<<new_v.capacity()
            //        <<" v1 size:"<<v1->size()<<" v2 size:"<<v2->size()<<std::endl;
            new_v.shrink_to_fit();
        }
        delete [] heap;
        delete [] lin;
        return new_v;  
    }
}
#endif