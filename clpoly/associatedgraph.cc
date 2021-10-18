/**
 * @file associatedgraph.cc
 * @author  李昊坤 (ker@pm.me)
 * @brief 对多项式集合生成associatedgraph, 以及求相关弦性操作等。
 * 
 */

#include<clpoly/associatedgraph.hh>

namespace clpoly{ 
    void __MCS(const std::vector<std::set<uint64_t>> &g,std::vector<uint64_t>& l)
    {
        std::vector<uint64_t> w,l_;
        w.resize(g.size());
        l.resize(g.size());
        l_.resize(g.size());
        uint64_t wmax,windex;
        bool tmp_b;
        for (uint64_t i=g.size();i>0;--i)
        {
            tmp_b=true;
            for (uint64_t j=0;j<w.size();++j)
            {
                if (!l_[j] && (tmp_b || w[j]>wmax))
                {
                    tmp_b=false;
                    wmax=w[j];
                    windex=j;
                }
                else
                {
                    if (!l_[j] && w[j]==wmax && std::rand()>RAND_MAX /2 )
                    {
                        windex=j;
                    }
                }
                
            }
            l[i-1]=windex;
            l_[windex]=i;
            for (auto & j :g[windex])
            {
                if (!l_[j])
                    ++w[j];
            }
        }
    }
    void __MCS_M(std::vector<std::set<uint64_t>> &g,std::vector<uint64_t>& l) 
    {
        std::vector<uint64_t> w,tmp_1,tmp_2;
        std::vector<uint64_t> l_;
        std::vector<bool> tmp_3,tmp_4;
        std::srand(std::time(nullptr));
        w.resize(g.size());
        l.resize(g.size());
        l_.resize(g.size());
        tmp_1.resize(g.size());
        tmp_2.resize(g.size());
        tmp_3.resize(g.size());
        tmp_4.resize(g.size());
        uint64_t wmax,windex;
        bool tmp_b;
        for (uint64_t i=g.size();i>0;--i)
        {
            tmp_b=true;
            for (uint64_t j=0;j<w.size();++j)
            {
                if (!l_[j] && (tmp_b || w[j]>wmax))
                {
                    tmp_b=false;
                    wmax=w[j];
                    windex=j;
                }
                else
                {
                    if (!l_[j] && w[j]==wmax && std::rand()>RAND_MAX /2 )
                    {
                        windex=j;
                    }
                }
                
            }
            l[i-1]=windex;
            l_[windex]=i;
            if (i>1)
            {
                for (uint64_t j=0;j<w.size();++j)
                {
                    tmp_1[j]=-1;
                    tmp_3[j]=false;
                    tmp_4[j]=false;
                }
                for (auto &j:g[windex])
                    tmp_4[j]=true;
                tmp_1[windex]=0;
                tmp_3[windex]=true;
                tmp_2[0]=windex;
                uint64_t tmp_begin=0,tmp_end=1,tmp_5;
                while(tmp_begin!=tmp_end)
                {
                    uint64_t tmp_index=tmp_2[tmp_begin++ % w.size()];
                    tmp_3[tmp_index]=false;
                    for(auto &j:g[tmp_index])
                    {
                        if (j!=tmp_index && !l_[j])
                        {
                            //tmp_4=std::max(tmp_1[tmp_index],w[j]);
                            if (tmp_1[tmp_index]<w[j])
                            {
                                tmp_5=w[j];
                                if (!tmp_4[j] && !l_[j])
                                {
                                    g[j].insert(windex);
                                    g[windex].insert(j);
                                    tmp_4[j]=true;
                                    //std::cout<<i<<":"<<windex<<" "<<j<<std::endl;
                                }
                                
                            }
                            else
                            {
                                tmp_5=tmp_1[tmp_index];
                            }
                            if (tmp_1[j]==-1 || tmp_5<tmp_1[j])
                            {
                                if (!tmp_3[j])
                                {
                                    tmp_2[tmp_end++ % w.size()]=j;
                                    tmp_3[j]=true;
                                }
                                tmp_1[j]=tmp_5;
                            }
                        }
                    }

                }
                for (auto &j:g[windex])
                {
                    ++w[j];
                }
            }
        }
    }
    void _connected_branch_(uint64_t node,uint64_t branch,
                           const std::vector<std::set<uint64_t>> & g,
                           std::vector<std::vector<uint64_t>> &l1,
                           std::vector<bool>&l2)
    {
        for (auto &j:g[node])
            if (l2[j])
            {
                l1[branch].push_back(j);
                l2[j]=false;
                _connected_branch_(j,branch,g,l1,l2);
            }
    }
    std::vector<std::vector<uint64_t>> _connected_branch(const std::vector<std::set<uint64_t>> & g)
    {
        std::vector<std::vector<uint64_t>> l1;
        std::vector<bool> l2;
        l2.resize(g.size());
        for (auto i=l2.begin();i!=l2.end();++i)
            *i=true;
        for (uint64_t i=0;i<g.size();++i)
            if (l2[i])
            {
                l1.push_back({});
                _connected_branch_(i,l1.size()-1,g,l1,l2);
            }
        return l1;
    }
}
