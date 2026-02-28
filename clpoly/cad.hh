/**
 * @file cad.hh
 * @author ntimesp(nxp@mail.ustc.edu.cn) 李昊坤(ker@pm.me)
 * @brief  定义CAD相关函数
 */
#ifndef CLPOLY_CAD_HH
#define CLPOLY_CAD_HH

#include <clpoly/resultant.hh>
#include <clpoly/polynomial_gcd.hh>
#include <clpoly/realroot.hh>
// todo: 设计好cell类后添加
// #include <clpoly/cell.hh>

namespace clpoly{
    enum class projection_method {
        MCCALLUM,    // McCallum 投影算子（默认）
        LAZARD      // Lazard 投影算子
    };

    // 计算多项式集合 F 关于 x 的 conts 集和 prims 集
    // 保证 x 是 F 中可能出现的变元里最小的
    // 只处理 first_var == x 的多项式
    // 常数(包括0)多项式会被自动跳过
    template<class var_order>
    std::pair<std::vector<polynomial_<ZZ,lex_<var_order>>>,
            std::vector<polynomial_<ZZ,lex_<var_order>>>>
    __conts_prims_polys_var(const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys,const variable &x)
    {
        std::vector<polynomial_<ZZ,lex_<var_order>>> conts;
        std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
        for (const auto& poly : polys) 
        {
            if (is_number(poly))
                continue;
            // first_var 必须大于 x
            // 常量的 first_var 为 variable() 顺序不确定，因此要先排除常量
            assert(get_first_var(poly) == x || poly.comp(x, get_first_var(poly)));
            if (get_first_var(poly) == x)
            {
                auto f_cont=cont(poly);
                // cont 不能是0
                assert(!f_cont.empty());
                prims.push_back(poly/f_cont);
                if (!is_number(f_cont))
                    conts.push_back(std::move(f_cont));
            }
            else
            {
                conts.push_back(poly);
            }
        }
        return {conts,prims};
    }


    // 字典序的投影算子
    // 需要传入x, 因为某些多项式可能不含x
    // 假设x是多项式包含的所有变量中var_order最小的变量
    template <class var_order>
    std::vector<polynomial_<ZZ,lex_<var_order>>> __project(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const clpoly::variable& x, 
        projection_method method = projection_method::LAZARD
    )
    {
        if (polys.empty())
            return {};
        if (method == projection_method::MCCALLUM) 
        {
            return __project_mccallum(polys, x);
        }
        else 
        {
            return __project_lazard(polys, x);
        }
    }

    // 字典序的投影算子
    // 需要传入x, 因为某些多项式可能不含x
    // 假设x是多项式包含的所有变量中 var_order 最小的变量
    // 新版本：先算 cont, prim 再算 squarefreebasis
    // todo: 用一个开关确定, polys是不是本身已经是 prims 并且 square free
    template <class var_order>
    std::vector<polynomial_<ZZ,lex_<var_order>>> __project_mccallum(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const clpoly::variable& x
    )
    {
        
        // 计算cont和prim
        // std::vector<polynomial_<ZZ,lex_<var_order>>> projs,prims_raw;
        auto [projs,prims_raw]=__conts_prims_polys_var(polys,x);

        // 用 squarefreebasis 计算无平方的基，调用 字典序 的实现
        auto [sqfree_prims, _] = squarefreebasis(prims_raw);

        // 注意sqfree_prims 可能 不含 x, 过滤掉不含 x 的多项式
        std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
        prims.reserve(sqfree_prims.size());

        for (auto& poly : sqfree_prims) {
            if (is_number(poly))
                continue;

            if (get_first_var(poly) == x) {
                prims.push_back(std::move(poly));
            } else {
                // 例如 (x-1)*g(y) 拆出来的 g(y)，只含低维变量，下传
                projs.push_back(std::move(poly));
            }
        }

        // 计算 coeff 和 discriminant, 保证每个多项式都含x
        for (const auto& poly : prims) {
            // 保证每个多项式都含x
            assert(get_first_var(poly)==x);

            // 调用 字典序 的 coeff
            auto coF = coeff(poly);
            for (auto& c : coF) {
                if (!is_number(c))
                    projs.push_back(std::move(c));
            }
            // discriminant 只有 一般序 和 univariate_priority_order 的实现，调用一般序的实现
            // todo: 以后可以添加 字典序 的 discriminant 实现
            auto disc = discriminant(poly, x);
            assert(!disc.empty());
            if (!is_number(disc))
                projs.push_back(std::move(disc));
        }
        //计算每对多项式的 resultant
        for (size_t i = 0; i < prims.size(); ++i) {
            for (size_t j = i + 1; j < prims.size(); ++j) {
                // resultant 只有 一般序 和 univariate_priority_order 的实现，调用一般序的实现
                // todo: 以后可以添加 字典序 的 resultant 实现
                auto res = resultant(prims[i], prims[j], x);
                assert(!res.empty());
                if (!is_number(res))
                    projs.push_back(std::move(res));
            }
        }
        return projs;
    }

    // 快速判断多项式的 first_var 是不是多项式的因子
    // 常量的 first_var 为 variable()，认为 variable() 不是常数的因子, 返回 false
    // 例子：假设变量序 x<y<z
    // 3 的 first_var 是 variable()，输出 false
    // x*y+y 的 first_var 是 x，输出 false  
    // x*y+x 的 first_var 是 x，输出 true
    template <class var_order>
    bool __has_factor_first_var(
        const polynomial_<ZZ,lex_<var_order>>& poly
    )
    {
        auto x = get_first_var(poly);
        // 只需要检查最后一项的 first_var 是否为x
        return !is_number(poly) && get_first_var(poly.back().first) == x;
    }
    
    // 计算末尾项系数，只有lazard投影算子会用到
    // 如果整除x，则 tailcoeff 为0
    // 常数的 tailcoeff 为 目前实现为常数本身
    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> __tailcoeff_lazard(const polynomial_<ZZ,lex_<var_order>> &F_)
    {
        polynomial_<ZZ,lex_<var_order>>  tc(F_.comp_ptr());
        auto v=get_first_var(F_);
        // 只需要从后向前把不含v的项加入tc，直到遇到含v的项为止
        for (auto it = F_.end(); it != F_.begin(); ) {
            --it;
            if (!it->first.empty() && it->first.front().first == v) {
                break;
            }
            tc.push_back(*it);
        }
        // 需要把tc反转回来
        std::reverse(tc.begin(), tc.end());
        return tc;
    }

    // 字典序的投影算子，只能被CAD调用
    // 需要传入x, 因为某些多项式可能不含x
    // 假设x是多项式包含的所有变量中 var_order 最小的变量
    // 新版本：先算 cont, prim 再算 squarefreebasis
    // todo: 用一个开关确定, polys是不是本身已经是 prims 并且 square free
    template <class var_order>
    std::vector<polynomial_<ZZ,lex_<var_order>>> __project_lazard(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const clpoly::variable& x
    )
    {
        // 计算cont和prim
        // std::vector<polynomial_<ZZ,lex_<var_order>>> projs,prims_raw;
        auto [projs,prims_raw]=__conts_prims_polys_var(polys,x);

        // 用 squarefreebasis 计算无平方的基，调用 字典序 的实现
        auto [sqfree_prims, _] = squarefreebasis(prims_raw);

        // 注意sqfree_prims 可能 不含 x, 过滤掉不含 x 的多项式
        // 如果有因子x, 需要分解   
        std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
        prims.reserve(sqfree_prims.size());

        for (auto& poly : sqfree_prims) {
            if (is_number(poly))
                continue;

            if (get_first_var(poly) == x) {
                // 如果有因子x, 需要分解
                if (__has_factor_first_var(poly))
                {
                    // 强制类型转换 x
                    polynomial_<ZZ,lex_<var_order>> x_factor(poly.comp_ptr());
                    x_factor.push_back({{{x,1}},1});
                    prims.push_back(std::move(x_factor));

                    // 原地做除法 poly /= x, 把首变量的次数减1, 如果次数为1则删除该变量
                    for (auto& term : poly) {
                        if (term.first.front().second > 1) {
                            // 下面是一段危险操作！
                            term.first.front().second -= 1;
                            term.first.deg() -= 1;   // 同时维护总次数
                        } 
                        else 
                        {
                            // erase 里已维护 __deg
                            term.first.erase(term.first.begin());
                        }
                        // 确保修改没有问题
                        assert(term.first.is_normal());
                    }
                }
                if (!is_number(poly))
                    prims.push_back(std::move(poly));
            } 
            else 
            {
                // 例如 (x-1)*g(y) 拆出来的 g(y)，只含低维变量，下传
                projs.push_back(std::move(poly));
            }
        }

        // 计算 leadcoeff, tailcoeff 和 discriminant, 保证每个多项式都含x
        for (const auto& poly : prims) {
            // 保证每个多项式都含x
            assert(get_first_var(poly)==x);

            // 调用 字典序 的 leadcoeff
            auto lc = leadcoeff(poly);
            if (!is_number(lc))
                projs.push_back(std::move(lc));
            // 调用 字典序 的 tailcoeff
            auto tc = __tailcoeff_lazard(poly);
            if (!is_number(tc))
                projs.push_back(std::move(tc));
            // discriminant 只有 一般序 和 univariate_priority_order 的实现，调用一般序的实现
            // todo: 以后可以添加 字典序 的 discriminant 实现
            auto disc = discriminant(poly, x);
            assert(!disc.empty());
            if (!is_number(disc))
                projs.push_back(std::move(disc));
        }

        //计算每对多项式的 resultant
        for (size_t i = 0; i < prims.size(); ++i) {
            for (size_t j = i + 1; j < prims.size(); ++j) {
                // resultant 只有 一般序 和 univariate_priority_order 的实现，调用一般序的实现
                // todo: 以后可以添加 字典序 的 resultant 实现
                auto res = resultant(prims[i], prims[j], x);
                assert(!res.empty());
                if (!is_number(res))
                    projs.push_back(std::move(res));
            }
        }
        return projs;
    }

    // 完整的投影, 目前不走分支
    // 保证 vars 符合多项式的 var_order
    // 计算第i层的投影的步骤
    // step1 计算 cont 和 prim
    // step2 计算 prim 的 squarefreebasis, 过滤掉不含 x_i 的多项式, 作为第l层的投影多项式
    // step3 只把过滤后的 squarefreebasis 传给 __project, 计算单次投影
    // 输出: 每层的投影多项式，确保第i层所有多项式都包含vars[i]
    template <class var_order>
    std::vector<std::vector<polynomial_<ZZ,lex_<var_order>>>> __project_full(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const std::vector<clpoly::variable>& vars, 
        projection_method method = projection_method::LAZARD
    )
    {
        // 没有变量，说明全是常数，返回空
        if(vars.empty()) return {};

        // 保存计算结果
        std::vector<std::vector<polynomial_<ZZ,lex_<var_order>>>> allprojs;
        // current_set 表示用于本层计算 cont/prim 的多项式集合（初始为输入 polys）
        std::vector<polynomial_<ZZ,lex_<var_order>>> current_set = polys;
        // 只需要投影前 n-1 个变量
        for(size_t i=0;i<vars.size();i++)
        {
            const auto &x = vars[i];
            // step1 计算 cont 和 prim
            auto tmp=__conts_prims_polys_var(current_set,x);
            auto conts = std::move(tmp.first);
            auto prims_raw = std::move(tmp.second);
            // step2: 用 squarefreebasis 计算无平方的基，调用 字典序 的实现
            auto [sqfree_prims, _] = squarefreebasis(prims_raw);
            // 过滤掉不含 x_l 的多项式, 作为第l层的投影多项式
            std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
            prims.reserve(sqfree_prims.size());
            conts.reserve(conts.size()+sqfree_prims.size());

            for (auto& poly : sqfree_prims) {
                if (is_number(poly))
                    continue;

                if (get_first_var(poly) == x) {
                    prims.push_back(std::move(poly));
                } 
                else 
                {
                    // 例如 (x-1)*g(y) 拆出来的 g(y)，只含低维变量，下传
                    conts.push_back(std::move(poly));
                }
            }
            // 把该层的 prims 作为第 l 层的投影多项式（所有元素都含 x）
            allprojs.push_back(std::move(prims));

            // 如果这是最后一层，不用再投影下一层，循环结束
            if (i + 1 >= vars.size()) {
                break;
            }

            // 计算投影
            current_set=__project(allprojs.back(), x, method);

            // 合并投影和conts, 作为下一层的current_set
            current_set.reserve(current_set.size() + conts.size());
            current_set.insert(current_set.end(),
                                std::make_move_iterator(conts.begin()),
                                std::make_move_iterator(conts.end()));
        }
        
        return allprojs; 
    }

    // 计算一维多项式没两个根之间的一个有理数
    template <class comp>
    std::vector<clpoly::QQ> sample_open_intervals(
        const std::vector<polynomial_<ZZ,comp>>& polys
    )
    {
        // 计算开区间样本点
        std::vector<clpoly::QQ> sample_points;

        // 调用 一般序 接口，内部会先全转成单变元多项式类型
        auto [roots,_] = realroot(polys);
        // 可能没有根
        if (roots.empty()) return {};
        
        auto root=roots.begin();
        sample_points.push_back(root->left() - 1);
        auto preright=root->right();
        root++;
        for (; root != roots.end(); ++root) {
            assert(preright < root->left());
            sample_points.push_back((preright + root->left()) / 2);
            preright = root->right();
        }
        sample_points.push_back(preright + 1);
        
        return sample_points;    
    }

    // // 字典序的 open cad
    // // 保证 vars 符合多项式的 var_order
    // // 输入: 多项式集, 变量序
    // // 输出: 全部开胞腔，每个胞腔中一个有理样本点
    // // 然后提升
    // template <class var_order>
    // std::pair<std::vector<cell<ZZ,lex_<var_order>>>,std::vector<std::vector<QQ>>>
    // __open_cad(
    //     const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
    //     const std::vector<variable>& vars
    // )
    // {
    //     // 先计算完整的投影, 使用默认的lazard投影算子
    //     auto allprojs=__project_full(polys,vars);

    //     // 提升

    // }

    //一般序，先转化为字典序 lex_<custom_var_order>，其中custom_var_order md({x})
    // 然后调用字典序方法
    template <class comp>
    std::vector<polynomial_<ZZ,comp>> project(
        const std::vector<polynomial_<ZZ,comp>>& polys, 
        const clpoly::variable& x, 
        projection_method method = projection_method::LAZARD
    )
    {
        std::vector<polynomial_<ZZ,comp>> projs;
        if (polys.empty())
            return projs;
        // 先转化为字典序
        std::vector<polynomial_<ZZ,lex_<custom_var_order>>> polys_lex;
        // custom_var_order 中 x 是最小变量
        lex_<custom_var_order> md({x});
        polynomial_<ZZ,lex_<custom_var_order>> p(&md);
        polys_lex.reserve(polys.size());
        for (auto &i:polys)
        {
            poly_convert(i,p);
            polys_lex.push_back(std::move(p));
        }
        // 计算字典序的投影
        auto projs_lex=__project(polys_lex,x,method);
        // 转化回原序
        projs.reserve(projs_lex.size());
        polynomial_<ZZ,comp> p1(polys.front().comp_ptr());
        for (auto &i:projs_lex)
        {
            poly_convert(i,p1);
            projs.push_back(std::move(p1));
        }
        return projs;
    }

}

#endif