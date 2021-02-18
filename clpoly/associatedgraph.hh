/**
 * @file associatedgraph.hh
 * @author  李昊坤 (ker@pm.me)
 * @brief 对多项式集合生成associatedgraph, 以及求相关弦性操作等。
 * 
 */

#ifndef CLPOLY_ASSOCIATEDGRAPH_HH
#define CLPOLY_ASSOCIATEDGRAPH_HH
#include <clpoly/polynomial.hh>
#include <vector>
#include <random>
#include <list>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>
#include <cassert>
namespace clpoly{ 
    /**
     * graph<node_type> 类是一个稀疏图类 
     * 内部有基本储存结构是将node_type映射到 index（uint64_t），以及一个基于index的连接表
     * index是从0到node的size-1
     * 构造方法只有 默认构造
     * 通过 add_node 增加节点
     * 通过 add_edge 增加边
     * nodes  mnode  adjacency_list 可以访问基本结构但不能修改
     * node 可以获取对应index的node
     * index 可以获取对应node的index
     * neighborhood_index 可以获取 node 相邻的 node 的 index
     * is_edge 可以查看边是否存在
     * size 可以查看图节点大小,同时也是图节点index上界
     **/
    template<class node_type>
    class graph
    {
        private:
            std::vector<node_type> __nodes;
            std::map<node_type,uint64_t> __mnode;
            std::vector<std::set<uint64_t>> __adjacency_list;  
        public:
            typedef node_type nodetype;

            
            constexpr const std::vector<node_type> & nodes() const {return this->__nodes;}
            constexpr const std::map<node_type,uint64_t> & mnode() const {return this->__mnode;}
            constexpr const std::vector<std::set<uint64_t>> & adjacency_list() const {return this->__adjacency_list;}  
            constexpr const node_type & node(uint64_t i)const {return this->__nodes[i];}
            constexpr const std::set<uint64_t> & neighborhood_index(const node_type & v) const {
                auto search=this->__mnode.find(v);
                assert (search!=this->__mnode.end());
                return this->__adjacency_list[search->second];
            }
            constexpr uint64_t index(const node_type & v) const {
                auto search=this->__mnode.find(v);
                assert (search!=this->__mnode.end());
                return search->second; 
            }

            constexpr auto size() const {return this->__nodes.size();}
            
            inline void add_node(node_type c)
            {
                auto search=this->__mnode.find(c);
                if (search==this->__mnode.end())
                {
                    this->__nodes.push_back(std::move(c));
                    this->__mnode[c]=this->__nodes.size()-1;
                    this->__adjacency_list.push_back({});
                }
            }

            void add_edge(const node_type &a,const node_type &b,bool directed=false)
            {
                uint64_t na,nb;
                auto search=this->__mnode.find(a);
                if (search==this->__mnode.end())
                {
                    this->add_node(a);
                    na=this->size()-1;
                }
                else{
                    na=search->second;
                }
                search=this->__mnode.find(b);
                if (search==this->__mnode.end())
                {
                    this->add_node(b);
                    nb=this->size()-1;
                }
                else{
                    nb=search->second;
                }
                this->__adjacency_list[na].insert(nb);
                if (!directed)
                {
                    this->__adjacency_list[nb].insert(na);
                }
            }
            void add_edge_index(uint64_t na,uint64_t nb,bool directed=false)
            {
                this->add_edge(na,nb,directed);
            }

            void add_edge(uint64_t na,uint64_t nb,bool directed=false)
            {
                assert(na<__nodes.size() && nb<__nodes.size());
                this->__adjacency_list[na].insert(nb);
                if (!directed)
                {
                    this->__adjacency_list[nb].insert(na);
                }
            }

            bool is_edge(uint64_t na,uint64_t nb) const
            {
                assert(na<__nodes.size() && nb<__nodes.size());
                auto search1=this->__adjacency_list[na].find(nb);
                if (search1==this->__adjacency_list[na].end())
                    return false;
                return true;
            }
            bool is_edge(const node_type &a,const node_type &b) const
            {
                uint64_t na,nb;
                auto search=this->__mnode.find(a);
                if (search==this->__mnode.end())
                {
                    return false;
                }
                else{
                    na=search->second;
                }
                search=this->__mnode.find(b);
                if (search==this->__mnode.end())
                {
                    return false;
                }
                else{
                    nb=search->second;
                }
                auto search1=this->__adjacency_list[na].find(nb);
                if (search1==this->__adjacency_list[na].end())
                    return false;
                return true;
            }
            
            template<class T>
            friend  void chordal_completion(graph<T>& gout,const graph<T>& gin,std::vector<uint64_t>& l);
            
    };

    template<class node_type>
    std::ostream& operator<<  (std::ostream& stream, const graph<node_type>& G) {
        stream<<"{\n";
        for(auto &i:G.mnode())
        {
            stream<<" "<<i.first<<":";
            const std::set<uint64_t> & S=G.adjacency_list()[i.second];
            for (auto & j:S)
            {
                stream<<" "<<G.nodes()[j];
            }
            stream<<"\n";
        }        
        stream<<"}";
        return stream;
    }

    template<class Tc,class comp>
    graph<variable> associatedgraph(const std::vector<polynomial_<Tc,comp>> & polys)
    {
        graph<variable> G;
        for(auto & p:polys)
        {
            auto vs=p.variables();
            for (auto i=vs.begin();i!=vs.end();++i)
            {
                G.add_node(i->first);
                for (auto j=vs.begin();j!=i;++j)
                    G.add_edge(i->first,j->first);
            }
        }
        return G;
    }
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

    
    template<class node_type>
    void chordal_completion(graph<node_type>& gout,const graph<node_type>& gin,std::vector<uint64_t>& l)
    {
        gout=gin;
        __MCS_M(gout.__adjacency_list,l);
    }
    template<class node_type>
    std::vector<node_type> chordal_completion(graph<node_type>& gout,const graph<node_type>& gin)
    {
        gout=gin;
        std::vector<uint64_t> l1;
        chordal_completion(gout,gin,l1);
        std::vector<node_type> l;
        l.resize(l1.size());
        for(auto i=0;i<l.size();i++)
            l[i]=gin.node(l1[i]);
        return l;
    }
    template<class node_type>
    graph<node_type> chordal_completion(const graph<node_type>& g)
    {
        graph<node_type>  gout=g;
        std::vector<uint64_t> l;
        chordal_completion(gout,g,l);
        return gout;
    }
    template<class node_type>
    std::vector<node_type> perfect_elimination_ordering(const graph<node_type>& gin)
    {
        std::vector<uint64_t> l1;
        __MCS(gin.adjacency_list(),l1);
        std::vector<node_type> l;
        l.resize(l1.size());
        for(auto i=0;i<l.size();i++)
            l[i]=gin.node(l1[i]);
        return l;
    }
    template<class node_type>
    double graph_diff_score(const graph<node_type>& g1,const graph<node_type>& g2)
    {
        double cross=0;
        // for (auto &i:g1.adjacency_list())
        //     sum+=i.size();
        // for (auto &i:g2.adjacency_list())
        //     sum+=i.size();
        for(uint64_t i=0;i<g1.size();++i)
            for(auto &j:g1.adjacency_list()[i])
            {
                if (g2.is_edge(g1.node(i),g1.node(j)))
                    ++cross;
            }
        std::set<node_type> nodes;
        for (auto &i:g1.nodes())
            nodes.insert(i);
        for (auto &j:g1.nodes())
            nodes.insert(j);
        std::set<node_type> neighbor;
        double union_=0;
        for (auto &i:nodes)
        {
            neighbor.clear();
            auto search=g1.mnode().find(i);
            if (search!=g1.mnode().end())
            {
                for (auto &j:g1.adjacency_list()[search->second])
                    neighbor.insert(g1.node(j));   
                
            }
            search=g2.mnode().find(i);
            if (search!=g2.mnode().end())
            {
                for (auto &j:g2.adjacency_list()[search->second])
                    neighbor.insert(g2.node(j));                   
            }
            union_+=neighbor.size();
        }
        // std::cout<<cross<<" "<<union_<<std::endl;
        return 1-cross/union_;
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

    template<class node_type>
    std::vector<std::vector<node_type>> connected_branch(const graph<node_type> & g)
    {
        auto l=_connected_branch( g.adjacency_list());
        std::vector<std::vector<node_type>> l1;
        std::vector<node_type> l2;
        l1.reserve(l.size());
        for (auto &i:l)
        {
            l2.reserve(i.size());
            for (auto &j:i)
                l2.push_back(g.node(j));
            l1.push_back(std::move(l2));
        }
        return l1;
    }
    
    template<class node_type>
    graph<node_type> connected_branch_graph(const graph<node_type> & gin)
    {
        graph<node_type> gout=gin;
        auto l=_connected_branch(gout.adjacency_list());
        for (auto &l1:l)
            for (auto i=l1.begin();i!=l1.end();++i)
            {
                auto j=i;++j;
                for (;j!=l1.end();++j)
                    gout.add_edge_index(*i,*j);
            }
        return gout;
    }



    template<class node_type>
    graph<node_type> elimination_game(const graph<node_type>& g,const std::vector<node_type> & l)
    {
        graph<node_type>  gout=g;
        std::vector<bool> h;
        h.resize(g.size());
        for(auto i=h.begin();i!=h.end();++i)
            *i=true;
        for (auto &i:l)
        {
            auto i1=gout.index(i);
            h[i1]=false;
            auto & l1=gout.adjacency_list()[i1];
            for (auto j1=l1.begin();j1!=l1.end();++j1)
                if (h[*j1])
                {
                    auto j2=j1;
                    ++j2;
                    for (;j2!=l1.end();++j2)
                    {
                        if (h[*j2])
                        {
                            gout.add_edge_index(*j1,*j2);
                        }
                    }
                }
        }
        return gout;
    }
    template<class node_type>
    uint64_t elimination_height(const graph<node_type>& g,const std::vector<node_type> & l)
    {
        graph<node_type>  gout=g;
        std::vector<bool> h;
        h.resize(g.size());
        std::vector<uint64_t> height;
        std::map<uint64_t,uint64_t> lm;
        height.resize(g.size());
        uint64_t treeheight=0;
        for (uint64_t i=0;i<g.size();++i)
        {
            lm[gout.index(l[i])]=i;
            height[i]=0;
            h[i]=true;
        }
        for (auto &i:l)
        {
            
            auto i1=gout.index(i);
            h[i1]=false;
            if (height[i1]>treeheight)
                treeheight=height[i1];
            auto & l1=gout.adjacency_list()[i1];
            uint64_t minc=g.size();
            uint64_t minc_h=g.size()+1;
            for (auto j1=l1.begin();j1!=l1.end();++j1)
                if (h[*j1])
                {
                    if (lm[*j1]<minc_h)
                    {
                        minc=*j1;
                        minc_h=lm[*j1];
                    }
                    auto j2=j1;
                    ++j2;
                    for (;j2!=l1.end();++j2)
                    {
                        if (h[*j2])
                        {
                            gout.add_edge_index(*j1,*j2);
                        }
                    }
                }
            if (minc!=g.size())
                if (height[i1]>=height[minc])
                    height[minc]=height[i1]+1;
                
        }
        return treeheight;
    }

    template<class Tc>
    std::pair<uint64_t,uint64_t> __polynomial_m_d(const std::vector<polynomial_<Tc>>& P,const std::vector<variable> & l)
    {
        std::vector<uint64_t> md;
        md.resize(l.size());
        for (auto &i:md)
            i=0;
        std::map<variable,uint64_t> ml;
        for (uint64_t i=0;i<l.size();++i)
            ml[l[i]]=i;
        uint64_t d=0;
        for(auto & p:P)
        {
            auto vs=p.variables();
            uint64_t m=l.size();
            for (auto i=vs.begin();i!=vs.end();++i)
            {
                if (ml[i->first]<m)
                    m=ml[i->first];
                if (i->second>d)
                    d= i->second;
            }
            ++md[m];

        }
        uint64_t m=0;
        for (auto &i:md)
            if (i>m)
                m=i;
        return {m,d};
    }

}
#endif