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
namespace clpoly{
    template<class T>
    inline bool zore_check(const T & c)
    {
        return c==0;
    }
    template<class T>
    inline void add_assignment(T & x,const T & y)
    {
        x+=y;
    } 
    template<class T>
    inline void sub_assignment(T & x,const T & y)
    {
        x-=y;
    } 
    template<class T>
    inline T negate (const T & m)
    {
        return -m;
    }
    template<class T>
    inline bool greater(const T & m1,const T & m2)
    {
        return m1>m2;
    }
    template<class T>
    inline bool less(const T & m1,const T & m2)
    {
        return m1<m2;
    }
    template<class T>
    inline bool  equal_to(const T & m1,const T & m2)
    {
        return m1==m2;
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
            return negate(v2);
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
}
#endif