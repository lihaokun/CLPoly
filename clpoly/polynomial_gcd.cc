/**
 * @file polynomial_gcd.cc
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义polynomial_gcd
*/
#include<clpoly/polynomial_gcd.hh>
namespace clpoly{

    // ---- GCDHEU: 启发式 GCD (Parisse 2002 Theorem 1) ----
    // 返回 true = 成功（result 为 primitive GCD），false = 失败（回退到模算法）
    static bool __gcdheu(
        upolynomial_<ZZ>& result,
        const upolynomial_<ZZ>& F,
        const upolynomial_<ZZ>& G)
    {
        // F, G 必须已 primitive、非空、size > 1

        // Step 1: 计算最大系数位长
        size_t bits_f = 0, bits_g = 0;
        for (auto& term : F) {
            size_t b = term.second.sizeinbase(2);
            if (b > bits_f) bits_f = b;
        }
        for (auto& term : G) {
            size_t b = term.second.sizeinbase(2);
            if (b > bits_g) bits_g = b;
        }

        // 阈值检查：FLINT 用 2*FLINT_BITS = 128
        if (bits_f + bits_g >= 128) return false;

        // Step 2: pack_bits (Parisse bound 要求 +3，FLINT 用 +6 做安全裕量)
        size_t mn = std::min(bits_f, bits_g);
        size_t mx = std::max(bits_f, bits_g);
        size_t pack_bits = std::max(mn + 6, mx + 1);

        // Step 3: 稠密化 — 将 F, G 转为 dense 系数数组 [c_0, c_1, ..., c_deg]
        int64_t deg_f = get_deg(F);
        int64_t deg_g = get_deg(G);

        std::vector<ZZ> fc(deg_f + 1), gc(deg_g + 1);
        for (auto& term : F) fc[term.first.deg()] = term.second;
        for (auto& term : G) gc[term.first.deg()] = term.second;

        // Step 4: Horner 求值 F(ξ), G(ξ)  其中 ξ = 2^pack_bits
        // 乘以 ξ 即左移 pack_bits 位，加系数直接 mpz_add（允许负数）
        mpz_t f_val, g_val, g_int;
        mpz_inits(f_val, g_val, g_int, NULL);

        // Pack F: Horner: f_val = (...((c[deg]) * ξ + c[deg-1]) * ξ + ...) * ξ + c[0]
        mpz_set_ui(f_val, 0);
        for (int64_t i = deg_f; i >= 0; --i) {
            mpz_mul_2exp(f_val, f_val, pack_bits);       // f_val <<= pack_bits (×ξ)
            if (fc[i].is_small()) {
                int64_t v = fc[i].get_val();
                if (v >= 0) mpz_add_ui(f_val, f_val, static_cast<uint64_t>(v));
                else        mpz_sub_ui(f_val, f_val, static_cast<uint64_t>(-(v + 1)) + 1);
            } else {
                mpz_add(f_val, f_val, fc[i].get_mpz_v());
            }
        }

        // Pack G
        mpz_set_ui(g_val, 0);
        for (int64_t i = deg_g; i >= 0; --i) {
            mpz_mul_2exp(g_val, g_val, pack_bits);
            if (gc[i].is_small()) {
                int64_t v = gc[i].get_val();
                if (v >= 0) mpz_add_ui(g_val, g_val, static_cast<uint64_t>(v));
                else        mpz_sub_ui(g_val, g_val, static_cast<uint64_t>(-(v + 1)) + 1);
            } else {
                mpz_add(g_val, g_val, gc[i].get_mpz_v());
            }
        }

        // Step 5: 整数 GCD（取绝对值确保非负）
        mpz_abs(f_val, f_val);
        mpz_abs(g_val, g_val);
        mpz_gcd(g_int, f_val, g_val);

        // Step 6: 对称 ξ-adic 重构 g_int → 多项式 H
        mpz_t mask, half_xi, coeff_mpz;
        mpz_inits(mask, half_xi, coeff_mpz, NULL);
        mpz_set_ui(mask, 1);
        mpz_mul_2exp(mask, mask, pack_bits);    // mask = ξ = 2^pack_bits
        mpz_set(half_xi, mask);
        mpz_tdiv_q_ui(half_xi, half_xi, 2);    // half_xi = ξ/2
        mpz_sub_ui(mask, mask, 1);             // mask = ξ - 1

        upolynomial_<ZZ> H;
        int64_t deg_bound = std::min(deg_f, deg_g);
        for (int64_t i = 0; i <= deg_bound && mpz_sgn(g_int) != 0; ++i) {
            mpz_and(coeff_mpz, g_int, mask);
            mpz_tdiv_q_2exp(g_int, g_int, pack_bits);
            // 对称表示: 若 r > ξ/2，则 r -= ξ 并向高位进位 +1
            if (mpz_cmp(coeff_mpz, half_xi) > 0) {
                mpz_sub(coeff_mpz, coeff_mpz, mask);   // r -= (ξ-1)
                mpz_sub_ui(coeff_mpz, coeff_mpz, 1);   // 合计 r -= ξ（mask = ξ-1）
                mpz_add_ui(g_int, g_int, 1);            // 进位
            }

            if (mpz_sgn(coeff_mpz) != 0)
                H.push_back({umonomial(i), ZZ(static_cast<mpz_srcptr>(coeff_mpz))});
        }

        mpz_clears(f_val, g_val, g_int, mask, half_xi, coeff_mpz, NULL);

        if (H.empty()) return false;

        // 降幂排列（upolynomial 约定高次在前）
        std::reverse(H.data().begin(), H.data().end());
        H.normalization();

        // Step 7: pp(H)
        ZZ h_cont = cont(H);
        if (h_cont != 1 && h_cont != -1)
            for (auto& term : H) term.second /= h_cont;
        if (!H.empty() && H.front().second < 0)
            for (auto& term : H) term.second = -term.second;

        // Step 8: 验证 H | F 且 H | G
        upolynomial_<ZZ> Q, R;
        pair_vec_div(Q.data(), R.data(), F.data(), H.data(), F.comp());
        if (!R.empty()) return false;

        pair_vec_div(Q.data(), R.data(), G.data(), H.data(), F.comp());
        if (!R.empty()) return false;

        result = std::move(H);
        return true;
    }

    upolynomial_<ZZ>  polynomial_GCD(upolynomial_<ZZ> G,upolynomial_<ZZ> F)
    {
         if (F.empty())
            return G;
        if (G.empty())
            return F;
        if (F.size()==1 || G.size()==1)
        {
            auto ptr=F.begin();
            auto m=ptr->first;
            auto c=ptr->second;
            for(++ptr;ptr!=F.end();++ptr)
            {
                m=gcd(m,ptr->first);
                c=gcd(c,ptr->second);
            }
            for(ptr=G.begin();ptr!=G.end();++ptr)
            {
                m=gcd(m,ptr->first);
                c=gcd(c,ptr->second);
            }
            upolynomial_<ZZ> Pout;
            Pout={{m,c}};
            return Pout;
            
        }

        ZZ f_cont=cont(F);
        ZZ g_cont=cont(G);
        ZZ cont_gcd=gcd(f_cont,g_cont);
        for (auto &i:F)
            i.second/=f_cont;
        for (auto &i:G)
            i.second/=g_cont;

        // GCDHEU 快速路径
        {
            upolynomial_<ZZ> heu_result;
            if (__gcdheu(heu_result, F, G)) {
                if (cont_gcd != 1) {
                    for (auto& term : heu_result) term.second *= cont_gcd;
                }
                return heu_result;
            }
        }

        uint64_t prime = UINT64_C(18446744073709551557);  // 2^64 - 59
        upolynomial_<ZZ> Pout_,tmp_Pout_,R;
        upolynomial_<Zp> Pout_mod,f_p,g_p;
        upolynomial_<ZZ> tmp;
        ZZ lc_gcd=gcd(F.begin()->second,G.begin()->second );
        Zp lc_gcd_p;
        ZZ Pout_prime;
        std::int64_t Pout_d=INT64_MAX;
        std::int64_t tmp_Pout_d=INT64_MAX;

        // CRT 循环：每轮用一个素数计算 modular GCD，然后合并。
        // 素数从 2^64-59 开始递减，所需素数个数由 Mignotte bound 决定：
        // O(d*h/64)（d=度数, h=系数 bit-length），远小于可用素数总数（≈2^58）。
        // 与 GCL Algorithm 7.1 的 while(true) 结构一致。
        while (1)
        {

            while (F.begin()->second % prime ==0 || G.begin()->second % prime ==0)
            {
                prime = prev_prime_64(prime);
            }
            f_p=polynomial_mod(F,prime);
            g_p=polynomial_mod(G,prime);
            lc_gcd_p=Zp(lc_gcd,prime);

            // std::cout<<"p:"<<prime<<std::endl;
            // std::cout<<"f_p:"<<f_p<<std::endl;
            // std::cout<<"g_p:"<<g_p<<std::endl;
            // std::cout<<"lc_gcd_p:"<<lc_gcd_p<<std::endl;

            tmp_Pout_d=__polynomial_GCD(Pout_mod,f_p,g_p,lc_gcd_p,Pout_d);
            if (tmp_Pout_d==-1)
            {
                prime = prev_prime_64(prime);
                continue;
            }
            
            // std::cout<<"poly_mod:"<<Pout_mod<<std::endl;
            
            if (tmp_Pout_d < Pout_d)
            {
                Pout_d=tmp_Pout_d;
                Pout_prime=prime;
                // poly_convert(Pout_mod,Pout_);
                Pout_.clear();
                Pout_.reserve(Pout_mod.size());
                for (auto & i:Pout_mod)
                    Pout_.push_back({i.first,i.second.number()});
                

                for (auto &i:Pout_)
                {
                    i.second%=Pout_prime;
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                }
                
                // std::cout<<"Pout_:"<<Pout_<<std::endl;
                
            }
            else
            {
                Zp tmp_inv(Pout_prime,prime);
                tmp_inv=tmp_inv.inv();
                tmp_Pout_.clear();
                tmp_Pout_.reserve(Pout_.size());
                auto Pout_ptr=Pout_.begin();
                auto Pout_end=Pout_.end();
                auto Pm_ptr=Pout_mod.begin();
                auto Pm_end=Pout_mod.end();
                while (Pout_ptr!=Pout_end && Pm_ptr!=Pm_end)
                {
                    if (F.comp(Pout_ptr->first,Pm_ptr->first))
                    {
                        tmp_Pout_.push_back({Pout_ptr->first,Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
                        ++Pout_ptr;
                    }
                    else
                    {
                        if (Pout_ptr->first==Pm_ptr->first)
                        {
                            tmp_Pout_.push_back({std::move(Pm_ptr->first),Pout_ptr->second+
                            (Pm_ptr->second.number()-Pout_ptr->second)*tmp_inv.number()*Pout_prime});
                            ++Pout_ptr;
                        }
                        else
                        {
                            tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                            
                        }
                        ++Pm_ptr;
                    }
                    
                }
                while (Pout_ptr!=Pout_end)
                {
                    tmp_Pout_.push_back({std::move(Pout_ptr->first),Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
                    ++Pout_ptr;
                }
                while (Pm_ptr!=Pm_end)
                {
                    tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                    ++Pm_ptr;
                }
                Pout_prime*=prime; 
                for (auto &i:tmp_Pout_)
                {
                    i.second%=Pout_prime;
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                    else if (i.second<-Pout_prime/2)
                    {
                        i.second+=Pout_prime;
                    }
                }
                
                // std::cout<<"Pout_:"<<tmp_Pout_<<std::endl;
                
                if (tmp_Pout_==Pout_)
                {
                    if (Pout_d>0)
                    {
                        auto cont_=cont(Pout_);
                        // std::cout<<"cont_:"<<cont_<<std::endl;
                        // pair_vec_div(tmp.data(),tmp_Pout_.data(),cont_.data(),F.comp());
                        tmp=tmp_Pout_/upolynomial_<ZZ>({{0,cont_}});
                        std::swap(tmp,tmp_Pout_);
                        pair_vec_div(tmp.data(),R.data(),F.data(),tmp_Pout_.data(),F.comp());
                        if (R.empty())
                        {
                            pair_vec_div(tmp.data(),R.data(),G.data(),tmp_Pout_.data(),F.comp());
                            if (R.empty())
                            {
                                // pair_vec_multiplies(tmp.data(),tmp_Pout_.data(),cont_gcd.data(),F.comp());
                                tmp=tmp_Pout_*upolynomial_<ZZ>({{0,cont_gcd}});
                                return tmp;
                            }
                        }
                    }
                    else
                    {
                        return {{0,cont_gcd}};
                    }
                }
                std::swap(tmp_Pout_.data(),Pout_.data());
                       
            }
            prime = prev_prime_64(prime);
        }
    }
}