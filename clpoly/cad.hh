/**
 * @file cad.hh
 * @author ntimesp(nxp@mail.ustc.edu.cn) 李昊坤(ker@pm.me)
 * @brief  定义CAD相关函数
 */
#ifndef CLPOLY_CAD_HH
#define CLPOLY_CAD_HH

#include <clpoly/resultant.hh>
#include <clpoly/polynomial_gcd.hh>
#include <clpoly/polynomial_factorize.hh>
#include <clpoly/realroot.hh>
#include <clpoly/cad_tree.hh>
#include <clpoly/polynomial_convert.hh>
#include <map>
#include <set>
#include <stdexcept>

namespace clpoly{
    enum class projection_method {
        MCCALLUM,    // McCallum 投影算子（默认）
        LAZARD       // Lazard 投影算子
    };

    // 投影多项式基的生成方法
    enum class basis_computation_method {
        SQUAREFREE,  // 无平方基（默认，目前更快）
        FACTOR       // 不可约因子基（更彻底，目前更慢）
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
        projection_method method = projection_method::LAZARD,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
    )
    {
        if (polys.empty())
            return {};
        if (method == projection_method::MCCALLUM) 
        {
            return __project_mccallum(polys, x, basis_method);
        }
        else 
        {
            return __project_lazard(polys, x, basis_method);
        }
    }

    // McCallum 投影算子核心实现
    // 输入 prims 必须已经过 basis 计算，且每个多项式都含 x
    // 结果直接追加到 projs
    template <class var_order>
    void __project_mccallum_prims_to(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& prims,
        const clpoly::variable& x,
        std::vector<polynomial_<ZZ,lex_<var_order>>>& projs  // 直接追加
    )
    {
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
    }

    // 字典序的投影算子
    // 需要传入x, 因为某些多项式可能不含x
    // 假设x是多项式包含的所有变量中 var_order 最小的变量
    // 新版本：先算 cont, prim 再算 basis
    template <class var_order>
    std::vector<polynomial_<ZZ,lex_<var_order>>> __project_mccallum(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const clpoly::variable& x,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
    )
    {
        // 计算cont和prim
        auto [projs, prims_raw] = __conts_prims_polys_var(polys, x);

        // 根据 basis_method 选择使用 squarefreebasis 或 factorbasis，调用 字典序 的实现
        std::vector<polynomial_<ZZ,lex_<var_order>>> basis_prims;
        if (basis_method == basis_computation_method::FACTOR) {
            // 用 factorbasis 计算不可约互素基
            auto [irr_prims, _] = factorbasis(prims_raw);
            basis_prims = std::move(irr_prims);
        } else {
            // 用 squarefreebasis 计算无平方的基（默认）
            auto [sqfree_prims, _] = squarefreebasis(prims_raw);
            basis_prims = std::move(sqfree_prims);
        }

        // 注意 basis_prims 可能 不含 x, 过滤掉不含 x 的多项式
        std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
        prims.reserve(basis_prims.size());

        for (auto& poly : basis_prims) {
            if (is_number(poly))
                continue;

            if (get_first_var(poly) == x) {
                prims.push_back(std::move(poly));
            } else {
                // 例如 (x-1)*g(y) 拆出来的 g(y)，只含低维变量，下传
                projs.push_back(std::move(poly));
            }
        }

        __project_mccallum_prims_to(prims, x, projs);
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

    // Lazard 投影算子核心实现
    // 输入 prims 必须已经过 basis 计算，且每个多项式都含 x
    // 前置条件：如果 basis 是 SQUAREFREE，调用方必须先处理 x 因子
    template <class var_order>
    void __project_lazard_prims_to(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& prims,
        const clpoly::variable& x,
        std::vector<polynomial_<ZZ,lex_<var_order>>>& projs
    )
    {
        // 计算 leadcoeff, tailcoeff 和 discriminant, 保证每个多项式都含x
        for (const auto& poly : prims) {
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
            auto disc = discriminant(poly, x);
            assert(!disc.empty());
            if (!is_number(disc))
                projs.push_back(std::move(disc));
        }

        //计算每对多项式的 resultant
        for (size_t i = 0; i < prims.size(); ++i) {
            for (size_t j = i + 1; j < prims.size(); ++j) {
                auto res = resultant(prims[i], prims[j], x);
                assert(!res.empty());
                if (!is_number(res))
                    projs.push_back(std::move(res));
            }
        }
    }

    // 字典序的投影算子，只能被CAD调用
    // 需要传入x, 因为某些多项式可能不含x
    // 假设x是多项式包含的所有变量中 var_order 最小的变量
    // 新版本：先算 cont, prim 再算 basis
    template <class var_order>
    std::vector<polynomial_<ZZ,lex_<var_order>>> __project_lazard(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const clpoly::variable& x,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
    )
    {
        // 计算cont和prim
        auto [projs, prims_raw] = __conts_prims_polys_var(polys, x);

        // 根据 basis_method 选择使用 squarefreebasis 或 factorbasis，调用 字典序 的实现
        std::vector<polynomial_<ZZ,lex_<var_order>>> basis_prims;
        if (basis_method == basis_computation_method::FACTOR) {
            // 用 factorbasis 计算不可约互素基
            auto [irr_prims, _] = factorbasis(prims_raw);
            basis_prims = std::move(irr_prims);
        } else {
            // 用 squarefreebasis 计算无平方的基（默认）
            auto [sqfree_prims, _] = squarefreebasis(prims_raw);
            basis_prims = std::move(sqfree_prims);
        }

        // 注意 basis_prims 可能 不含 x, 过滤掉不含 x 的多项式
        // 如果有因子x, 需要分解（SQUAREFREE 模式时）
        std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
        prims.reserve(basis_prims.size());

        for (auto& poly : basis_prims) {
            if (is_number(poly))
                continue;

            if (get_first_var(poly) == x) {
                // 当使用 SQUAREFREE 时，如果有因子x, 需要手动分解
                // 当使用 FACTOR 时，factorbasis 已经分解了 x 因子
                if (basis_method == basis_computation_method::SQUAREFREE && 
                    __has_factor_first_var(poly))
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
                        } else {
                            // erase 里已维护 __deg
                            term.first.erase(term.first.begin());
                        }
                        // 确保修改没有问题
                        assert(term.first.is_normal());
                    }
                }
                if (!is_number(poly))
                    prims.push_back(std::move(poly));
            } else {
                // 例如 (x-1)*g(y) 拆出来的 g(y)，只含低维变量，下传
                projs.push_back(std::move(poly));
            }
        }

        __project_lazard_prims_to(prims, x, projs);
        return projs;
    }

    // 完整的投影, 目前不走分支
    // 保证 vars 符合多项式的 var_order
    // 计算第i层的投影的步骤
    // step1 计算 cont 和 prim
    // step2 计算 prim 的 basis, 过滤掉不含 x_i 的多项式, 作为第l层的投影多项式
    // step3 只把过滤后的 basis 传给 __project, 计算单次投影
    // 输出: 每层的投影多项式，确保第i层所有多项式都包含vars[i]
    template <class var_order>
    std::vector<std::vector<polynomial_<ZZ,lex_<var_order>>>> __project_full(
        const std::vector<polynomial_<ZZ,lex_<var_order>>>& polys, 
        const std::vector<clpoly::variable>& vars, 
        projection_method method = projection_method::LAZARD,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
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
            // step2: 根据 basis_method 选择使用 squarefreebasis 或 factorbasis, 调用 字典序 的实现
            std::vector<polynomial_<ZZ,lex_<var_order>>> basis_prims;
            if (basis_method == basis_computation_method::FACTOR) {
                auto [irr_prims, _] = factorbasis(prims_raw);
                basis_prims = std::move(irr_prims);
            } else {
                auto [sqfree_prims, _] = squarefreebasis(prims_raw);
                basis_prims = std::move(sqfree_prims);
            }
            // 过滤掉不含 x_l 的多项式, 作为第l层的投影多项式
            // 同时处理 x 因子（Lazard + SQUAREFREE 模式时）
            std::vector<polynomial_<ZZ,lex_<var_order>>> prims;
            prims.reserve(basis_prims.size());
            conts.reserve(conts.size()+basis_prims.size());

            for (auto& poly : basis_prims) {
                if (is_number(poly))
                    continue;

                if (get_first_var(poly) == x) {
                    // Lazard + SQUAREFREE 模式：处理 x 因子
                    if (method == projection_method::LAZARD && 
                        basis_method == basis_computation_method::SQUAREFREE &&
                        __has_factor_first_var(poly))
                    {
                        // 强制类型转换 x
                        polynomial_<ZZ,lex_<var_order>> x_factor(poly.comp_ptr());
                        x_factor.push_back({{{x,1}},1});
                        prims.push_back(std::move(x_factor));

                        // 原地做除法 poly /= x
                        for (auto& term : poly) {
                            if (term.first.front().second > 1) {
                                term.first.front().second -= 1;
                                term.first.deg() -= 1;
                            } else {
                                term.first.erase(term.first.begin());
                            }
                            assert(term.first.is_normal());
                        }
                    }
                    if (!is_number(poly))
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

            // 计算投影：prims 已经是 basis 且只含 x，直接调用核心实现避免重复计算 basis
            // conts 包含不含 x 的多项式，直接追加到投影结果
            current_set = std::move(conts);
            if (method == projection_method::MCCALLUM) {
                __project_mccallum_prims_to(allprojs.back(), x, current_set);
            } else {
                __project_lazard_prims_to(allprojs.back(), x, current_set);
            }
        }
        
        return allprojs; 
    }

    //一般序，先转化为字典序 lex_<custom_var_order>，其中custom_var_order md({x})
    // 然后调用字典序方法
    template <class comp>
    std::vector<polynomial_<ZZ,comp>> project(
        const std::vector<polynomial_<ZZ,comp>>& polys, 
        const clpoly::variable& x, 
        projection_method method = projection_method::LAZARD,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
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
        auto projs_lex=__project(polys_lex,x,method,basis_method);
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

    // ========== Open CAD ==========

    // 一次性将 realroot 输出的全局根索引映射为 cad_root
    // poly_mult_info: realroot 返回的第二个元素，poly_mult_info[i] = {(poly_idx, multiplicity), ...}
    // num_polys: 该层的多项式个数
    inline std::vector<cad_root> __make_cad_roots(
        const std::vector<std::vector<std::pair<uint64_t,uint64_t>>>& poly_mult_info,
        size_t num_polys
    )
    {
        // poly_root_count[p] = 多项式 p 已出现的不同实根数（不计重数）
        std::vector<size_t> poly_root_count(num_polys, 0);
        std::vector<cad_root> result(poly_mult_info.size());

        for (size_t i = 0; i < poly_mult_info.size(); ++i)
        {
            assert(!poly_mult_info[i].empty());
            // 选第一个多项式作为 cad_root 的标识
            size_t pidx = poly_mult_info[i][0].first;
            assert(pidx < num_polys);
            size_t local_idx = poly_root_count[pidx];
            result[i] = cad_root(pidx, local_idx);

            // 对该根所属的所有多项式递增计数，保证后续根的局部索引正确
            for (auto& pr : poly_mult_info[i])
            {
                assert(pr.first < num_polys);
                poly_root_count[pr.first]++;
            }
        }
        return result;
    }

    // 在指定 level 的 parent 下根据实根创建所有 Sector 节点
    // Level 0 时 parent_idx = SIZE_MAX（无 parent）
    template <class var_order>
    void __lift_open_level(
        cad_tree<var_order>& tree,
        size_t level,
        size_t parent_idx,
        const std::vector<uroot>& roots,
        const std::vector<std::vector<std::pair<uint64_t,uint64_t>>>& poly_mult_info
    )
    {
        // 根数为 0: 整条实数线是一个 Sector (-∞, +∞)
        if (roots.empty())
        {
            tree.add_sector(level, parent_idx, cad_root::neginf(), cad_root::inf(), QQ(0));
            return;
        }

        // 预计算所有根对应的 cad_root
        std::vector<cad_root> cad_roots = __make_cad_roots(poly_mult_info, tree.level_polys(level).size());
        assert(cad_roots.size() == roots.size());

        // (-∞, roots[0])
        tree.add_sector(level, parent_idx, cad_root::neginf(), cad_roots[0],
                        roots[0].left() - 1);

        // (roots[i], roots[i+1]) for i = 0..m-2
        for (size_t i = 0; i + 1 < roots.size(); ++i)
        {
            assert(roots[i].right() < roots[i+1].left());  // 隔离区间不重叠
            tree.add_sector(level, parent_idx, cad_roots[i], cad_roots[i+1],
                            (roots[i].right() + roots[i+1].left()) / 2);
        }

        // (roots[m-1], +∞)
        tree.add_sector(level, parent_idx, cad_roots.back(), cad_root::inf(),
                        roots.back().right() + 1);
    }

    // 将 _level_polys[level] 代入采样点，得到低维 ZZ 多项式
    // sample_path[k] 对应 level k 的变量（提升序）
    template <class var_order>
    std::vector<polynomial_<ZZ, lex_<var_order>>> __eval_level_polys(
        const cad_tree<var_order>& tree,
        size_t level,
        const std::vector<QQ>& sample_path
    )
    {
        using poly_type = polynomial_<ZZ, lex_<var_order>>;
        assert(sample_path.size() == level);

        // 构造代入映射: sample_path[k] 对应 level k 的变量
        std::map<variable, QQ> subst_map;
        for (size_t k = 0; k < sample_path.size(); ++k)
            subst_map[tree.level_var(k)] = sample_path[k];

        std::vector<poly_type> result;
        result.reserve(tree.level_polys(level).size());

        for (const auto& poly : tree.level_polys(level))
        {
            // ZZ 多项式代入 QQ 值 → QQ 多项式
            auto p_qq = assign<QQ, ZZ, QQ>(poly, subst_map);

            // QQ 多项式转 ZZ 多项式（乘公分母，不影响实根）
            poly_type p_zz(poly.comp_ptr());
            poly_convert(p_qq, p_zz);

            // 检测 nullification: 多项式在采样点处恒为零
            // 目前不处理此退化情况，直接报错
            if (p_zz.empty())
                throw std::runtime_error(
                    "__eval_level_polys: nullification detected at level " + std::to_string(level)
                    + " (polynomial vanishes identically at sample point)");

            result.push_back(std::move(p_zz));
        }
        return result;
    }

    // 计算 Open CAD（仅 Sector 节点）
    // 要求输入多项式为 lex_<var_order> 字典序
    template <class var_order>
    cad_tree<var_order> __open_cad(
        const std::vector<polynomial_<ZZ, lex_<var_order>>>& polys,
        const std::vector<variable>& vars,
        projection_method method = projection_method::LAZARD,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
    )
    {
        assert(!vars.empty());

        // 1. 投影
        auto allprojs = __project_full(polys, vars, method, basis_method);

        // 2. 构造树（构造函数内部反转为提升序）
        cad_tree<var_order> tree(vars, std::move(allprojs));
        size_t n = tree.num_levels();

        // 3. Level 0: _level_polys[0] 为一元多项式，无需代入
        auto [roots_0, info_0] = realroot(tree.level_polys(0));
        __lift_open_level(tree, 0, SIZE_MAX, roots_0, info_0);

        // 4. Level 1 到 Level n-1: 逐层提升
        for (size_t level = 1; level < n; ++level)
        {
            for (size_t parent_idx = 0; parent_idx < tree.level_size(level - 1); ++parent_idx)
            {
                auto eval_polys = __eval_level_polys(tree, level, tree.get_sample_point(level - 1, parent_idx));
                auto [roots_k, info_k] = realroot(eval_polys);
                __lift_open_level(tree, level, parent_idx, roots_k, info_k);
            }
        }

        return tree;
    }

    // 一般序的 open_cad，先转为字典序 lex_<custom_var_order> 再调用 __open_cad
    // vars 按 lex 序排列（vars[0] 最小）
    // 边界情况：vars 不能为空；polys 可以为空（返回整条实数线的 1 个 Sector）
    template <class comp>
    cad_tree<custom_var_order> open_cad(
        const std::vector<polynomial_<ZZ, comp>>& polys,
        const std::vector<variable>& vars,
        projection_method method = projection_method::LAZARD,
        basis_computation_method basis_method = basis_computation_method::SQUAREFREE
    )
    {
        if (vars.empty())
            throw std::invalid_argument("open_cad: vars must not be empty");

        // 检查 polys 中的变量必须全部出现在 vars 中
        {
            std::set<variable> var_set(vars.begin(), vars.end());
            for (const auto& vp : get_variables(polys))
                if (var_set.find(vp.first) == var_set.end())
                    throw std::invalid_argument(
                        "open_cad: polynomial contains variable '"
                        + vp.first.name() + "' not in vars");
        }

        // 构造字典序并转换
        lex_<custom_var_order> md(vars);
        polynomial_<ZZ, lex_<custom_var_order>> p(&md);
        std::vector<polynomial_<ZZ, lex_<custom_var_order>>> polys_lex;
        polys_lex.reserve(polys.size());
        for (const auto& i : polys)
        {
            poly_convert(i, p);
            polys_lex.push_back(std::move(p));
        }
        return __open_cad(polys_lex, vars, method, basis_method);
    }

}

#endif