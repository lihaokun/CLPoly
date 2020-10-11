/*
Module Name:
    associatedgraph.hh
Abstract:
    对多项式集合生成associatedgraph, 以及求相关弦性操作等。
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_ASSOCIATEDGRAPH_HH
#define CLPOLY_ASSOCIATEDGRAPH_HH
#include <clpoly/polynomial.hh>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>
#include <cassert>
namespace clpoly{ 
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

            constexpr auto size() const {return this->__nodes.size();}
            
            inline void add_node(node_type c)
            {
                this->__nodes.push_back(std::move(c));
                this->__mnode[c]=this->__nodes.size()-1;
                this->__adjacency_list.push_back({});
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
                assert(na<__nodes.size() && nb<__nodes.size());
                this->__adjacency_list[na].insert(nb);
                if (!directed)
                {
                    this->__adjacency_list[nb].insert(na);
                }
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







    
    

}
#endif