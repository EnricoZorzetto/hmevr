// Generated by rstantools.  Do not edit by hand.

/*
    hbevr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    hbevr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with hbevr.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_gev_nobounds_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_gev_nobounds");
    reader.add_event(88, 86, "end", "model_gev_nobounds");
    return reader;
}
template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const T1__& mu,
             const T2__& k,
             const T3__& psi, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 7;
        local_scalar_t__ N(DUMMY_VAR__);
        (void) N;  // dummy to suppress unused var warning
        stan::math::initialize(N, DUMMY_VAR__);
        stan::math::fill(N, DUMMY_VAR__);
        stan::math::assign(N,rows(y));
        current_statement_begin__ = 8;
        local_scalar_t__ inv_k(DUMMY_VAR__);
        (void) inv_k;  // dummy to suppress unused var warning
        stan::math::initialize(inv_k, DUMMY_VAR__);
        stan::math::fill(inv_k, DUMMY_VAR__);
        stan::math::assign(inv_k,inv(k));
        current_statement_begin__ = 9;
        if (as_bool((primitive_value(logical_lt(k, 0)) && primitive_value(logical_gt((max(subtract(y, mu)) / psi), -(inv_k)))))) {
            current_statement_begin__ = 10;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "k<0 and max(y-mu)/psi > -1/k; found k, psi =";
            errmsg_stream__ << k;
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 11;
        if (as_bool(logical_lte(psi, 0))) {
            current_statement_begin__ = 12;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "psi<=0; found psi =";
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 13;
        if (as_bool(logical_gt(stan::math::fabs(k), 1e-15))) {
            current_statement_begin__ = 15;
            return stan::math::promote_scalar<fun_return_scalar_t__>((((-(N) * stan::math::log(psi)) - ((1 + inv_k) * sum(stan::math::log1p(multiply(subtract(y, mu), (k / psi)))))) - sum(stan::math::exp(multiply(-(inv_k), stan::math::log1p(multiply(subtract(y, mu), (k / psi))))))));
        } else {
            current_statement_begin__ = 19;
            return stan::math::promote_scalar<fun_return_scalar_t__>(((-(N) * stan::math::log(psi)) - sum(add(divide(subtract(y, mu), psi), stan::math::exp(divide(minus(subtract(y, mu)), psi))))));
        }
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const T1__& mu,
             const T2__& k,
             const T3__& psi, std::ostream* pstream__) {
    return gev_lpdf<false>(y,mu,k,psi, pstream__);
}
struct gev_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const T1__& mu,
             const T2__& k,
             const T3__& psi, std::ostream* pstream__) const {
        return gev_lpdf(y, mu, k, psi, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lcdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const T1__& mu,
             const T2__& k,
             const T3__& psi, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 23;
        local_scalar_t__ inv_k(DUMMY_VAR__);
        (void) inv_k;  // dummy to suppress unused var warning
        stan::math::initialize(inv_k, DUMMY_VAR__);
        stan::math::fill(inv_k, DUMMY_VAR__);
        stan::math::assign(inv_k,inv(k));
        current_statement_begin__ = 24;
        if (as_bool((primitive_value(logical_lt(k, 0)) && primitive_value(logical_gt((max(subtract(y, mu)) / psi), -(inv_k)))))) {
            current_statement_begin__ = 25;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "k<0 and max(y-mu)/psi > -1/k; found k, psi =";
            errmsg_stream__ << k;
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 26;
        if (as_bool(logical_lte(psi, 0))) {
            current_statement_begin__ = 27;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "psi<=0; found psi =";
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 28;
        if (as_bool(logical_gt(stan::math::fabs(k), 1e-15))) {
            current_statement_begin__ = 29;
            return stan::math::promote_scalar<fun_return_scalar_t__>(sum(minus(stan::math::exp(multiply(-(inv_k), stan::math::log(add(1, multiply(subtract(y, mu), (k / psi)))))))));
        } else {
            current_statement_begin__ = 31;
            return stan::math::promote_scalar<fun_return_scalar_t__>(stan::math::exp(sum(stan::math::log(minus(stan::math::exp(divide(minus(subtract(y, mu)), psi)))))));
        }
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct gev_lcdf_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const T1__& mu,
             const T2__& k,
             const T3__& psi, std::ostream* pstream__) const {
        return gev_lcdf(y, mu, k, psi, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lccdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
              const T1__& mu,
              const T2__& k,
              const T3__& psi, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 35;
        local_scalar_t__ inv_k(DUMMY_VAR__);
        (void) inv_k;  // dummy to suppress unused var warning
        stan::math::initialize(inv_k, DUMMY_VAR__);
        stan::math::fill(inv_k, DUMMY_VAR__);
        stan::math::assign(inv_k,inv(k));
        current_statement_begin__ = 36;
        if (as_bool((primitive_value(logical_lt(k, 0)) && primitive_value(logical_gt((max(subtract(y, mu)) / psi), -(inv_k)))))) {
            current_statement_begin__ = 37;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "k<0 and max(y-mu)/psi > -1/k; found k, psi =";
            errmsg_stream__ << k;
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 38;
        if (as_bool(logical_lte(psi, 0))) {
            current_statement_begin__ = 39;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "psi<=0; found psi =";
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 40;
        if (as_bool(logical_gt(stan::math::fabs(k), 1e-15))) {
            current_statement_begin__ = 41;
            return stan::math::promote_scalar<fun_return_scalar_t__>(sum(log1m_exp(minus(stan::math::exp(multiply(-(inv_k), stan::math::log(add(1, multiply(subtract(y, mu), (k / psi))))))))));
        } else {
            current_statement_begin__ = 43;
            return stan::math::promote_scalar<fun_return_scalar_t__>(sum(log1m_exp(minus(stan::math::exp(divide(minus(subtract(y, mu)), psi))))));
        }
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct gev_lccdf_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
              const T1__& mu,
              const T2__& k,
              const T3__& psi, std::ostream* pstream__) const {
        return gev_lccdf(y, mu, k, psi, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__, class RNG>
typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
gev_rng(const T0__& mu,
            const T1__& k,
            const T2__& psi, RNG& base_rng__, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        current_statement_begin__ = 46;
        if (as_bool(logical_lte(psi, 0))) {
            current_statement_begin__ = 47;
            std::stringstream errmsg_stream__;
            errmsg_stream__ << "[psi<=0; found psi =";
            errmsg_stream__ << psi;
            throw std::domain_error(errmsg_stream__.str());
        }
        current_statement_begin__ = 48;
        if (as_bool(logical_gt(stan::math::fabs(k), 1e-15))) {
            current_statement_begin__ = 49;
            return stan::math::promote_scalar<fun_return_scalar_t__>((mu + ((psi / k) * (-(1) + stan::math::exp((-(k) * stan::math::log(-(stan::math::log(uniform_rng(0, 1, base_rng__))))))))));
        } else {
            current_statement_begin__ = 51;
            return stan::math::promote_scalar<fun_return_scalar_t__>((mu - (psi * stan::math::log(-(stan::math::log(uniform_rng(0, 1, base_rng__)))))));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct gev_rng_functor__ {
    template <typename T0__, typename T1__, typename T2__, class RNG>
        typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
    operator()(const T0__& mu,
            const T1__& k,
            const T2__& psi, RNG& base_rng__, std::ostream* pstream__) const {
        return gev_rng(mu, k, psi, base_rng__, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_gev_nobounds : public prob_grad {
private:
        int N;
        int Mgen;
        vector_d y;
        std::vector<double> pr_psi;
        std::vector<double> pr_mu;
        std::vector<double> pr_k;
public:
    model_gev_nobounds(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_gev_nobounds(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_gev_nobounds_namespace::model_gev_nobounds";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 55;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 56;
            context__.validate_dims("data initialization", "Mgen", "int", context__.to_vec());
            Mgen = int(0);
            vals_i__ = context__.vals_i("Mgen");
            pos__ = 0;
            Mgen = vals_i__[pos__++];
            check_greater_or_equal(function__, "Mgen", Mgen, 0);
            current_statement_begin__ = 57;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(N));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 58;
            validate_non_negative_index("pr_psi", "2", 2);
            context__.validate_dims("data initialization", "pr_psi", "double", context__.to_vec(2));
            pr_psi = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("pr_psi");
            pos__ = 0;
            size_t pr_psi_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < pr_psi_k_0_max__; ++k_0__) {
                pr_psi[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 59;
            validate_non_negative_index("pr_mu", "2", 2);
            context__.validate_dims("data initialization", "pr_mu", "double", context__.to_vec(2));
            pr_mu = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("pr_mu");
            pos__ = 0;
            size_t pr_mu_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < pr_mu_k_0_max__; ++k_0__) {
                pr_mu[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 60;
            validate_non_negative_index("pr_k", "2", 2);
            context__.validate_dims("data initialization", "pr_k", "double", context__.to_vec(2));
            pr_k = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("pr_k");
            pos__ = 0;
            size_t pr_k_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < pr_k_k_0_max__; ++k_0__) {
                pr_k[k_0__] = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 63;
            num_params_r__ += 1;
            current_statement_begin__ = 64;
            num_params_r__ += 1;
            current_statement_begin__ = 65;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_gev_nobounds() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 63;
        if (!(context__.contains_r("k")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable k missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("k");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "k", "double", context__.to_vec());
        double k(0);
        k = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(k);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable k: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 64;
        if (!(context__.contains_r("psi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable psi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("psi");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "psi", "double", context__.to_vec());
        double psi(0);
        psi = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, psi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable psi: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 65;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "mu", "double", context__.to_vec());
        double mu(0);
        mu = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 63;
            local_scalar_t__ k;
            (void) k;  // dummy to suppress unused var warning
            if (jacobian__)
                k = in__.scalar_constrain(lp__);
            else
                k = in__.scalar_constrain();
            current_statement_begin__ = 64;
            local_scalar_t__ psi;
            (void) psi;  // dummy to suppress unused var warning
            if (jacobian__)
                psi = in__.scalar_lb_constrain(0, lp__);
            else
                psi = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 65;
            local_scalar_t__ mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.scalar_constrain(lp__);
            else
                mu = in__.scalar_constrain();
            // model body
            current_statement_begin__ = 68;
            lp_accum__.add(normal_log<propto__>(mu, get_base1(pr_mu, 1, "pr_mu", 1), get_base1(pr_mu, 2, "pr_mu", 1)));
            current_statement_begin__ = 69;
            lp_accum__.add(lognormal_log<propto__>(psi, get_base1(pr_psi, 1, "pr_psi", 1), get_base1(pr_psi, 2, "pr_psi", 1)));
            current_statement_begin__ = 70;
            lp_accum__.add(normal_log<propto__>(k, get_base1(pr_k, 1, "pr_k", 1), get_base1(pr_k, 2, "pr_k", 1)));
            current_statement_begin__ = 71;
            lp_accum__.add(gev_lpdf<propto__>(y, mu, k, psi, pstream__));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("k");
        names__.push_back("psi");
        names__.push_back("mu");
        names__.push_back("yrep");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(Mgen);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_gev_nobounds_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double k = in__.scalar_constrain();
        vars__.push_back(k);
        double psi = in__.scalar_lb_constrain(0);
        vars__.push_back(psi);
        double mu = in__.scalar_constrain();
        vars__.push_back(mu);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 74;
            validate_non_negative_index("yrep", "Mgen", Mgen);
            std::vector<double> yrep(Mgen, double(0));
            stan::math::initialize(yrep, DUMMY_VAR__);
            stan::math::fill(yrep, DUMMY_VAR__);
            current_statement_begin__ = 75;
            validate_non_negative_index("log_lik", "N", N);
            std::vector<double> log_lik(N, double(0));
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 76;
            for (int i = 1; i <= Mgen; ++i) {
                current_statement_begin__ = 77;
                stan::model::assign(yrep, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            gev_rng(mu, k, psi, base_rng__, pstream__), 
                            "assigning variable yrep");
            }
            current_statement_begin__ = 79;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 80;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            gev_lpdf(rep_vector(get_base1(y, i, "y", 1), 1), mu, k, psi, pstream__), 
                            "assigning variable log_lik");
            }
            // validate, write generated quantities
            current_statement_begin__ = 74;
            size_t yrep_k_0_max__ = Mgen;
            for (size_t k_0__ = 0; k_0__ < yrep_k_0_max__; ++k_0__) {
                vars__.push_back(yrep[k_0__]);
            }
            current_statement_begin__ = 75;
            size_t log_lik_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
                vars__.push_back(log_lik[k_0__]);
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_gev_nobounds";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "k";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "psi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t yrep_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < yrep_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "yrep" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_k_0_max__ = N;
        for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "k";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "psi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "mu";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t yrep_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < yrep_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "yrep" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_k_0_max__ = N;
        for (size_t k_0__ = 0; k_0__ < log_lik_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_gev_nobounds_namespace::model_gev_nobounds stan_model;
#endif