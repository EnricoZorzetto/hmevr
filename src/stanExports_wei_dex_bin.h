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
namespace model_wei_dex_bin_namespace {
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
    reader.add_event(0, 0, "start", "model_wei_dex_bin");
    reader.add_event(111, 109, "end", "model_wei_dex_bin");
    return reader;
}
template <bool propto, typename T0__, typename T1__, typename T2__>
typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
exploc_lpdf(const T0__& x,
                const T1__& m,
                const T2__& s, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 16;
        local_scalar_t__ xmm(DUMMY_VAR__);
        (void) xmm;  // dummy to suppress unused var warning
        stan::math::initialize(xmm, DUMMY_VAR__);
        stan::math::fill(xmm, DUMMY_VAR__);
        current_statement_begin__ = 17;
        stan::math::assign(xmm, (x - m));
        current_statement_begin__ = 19;
        return stan::math::promote_scalar<fun_return_scalar_t__>(exponential_log(xmm, inv(s)));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__>
typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
exploc_lpdf(const T0__& x,
                const T1__& m,
                const T2__& s, std::ostream* pstream__) {
    return exploc_lpdf<false>(x,m,s, pstream__);
}
struct exploc_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
    operator()(const T0__& x,
                const T1__& m,
                const T2__& s, std::ostream* pstream__) const {
        return exploc_lpdf(x, m, s, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_wei_dex_bin : public prob_grad {
private:
        int M;
        int Mgen;
        int Nt;
        std::vector<std::vector<double> > y;
        std::vector<int> N;
        std::vector<double> cexploc0prior;
        std::vector<double> wexploc0prior;
        std::vector<double> cexpsc0prior;
        std::vector<double> wexpsc0prior;
        std::vector<double> pn0prior;
public:
    model_wei_dex_bin(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_wei_dex_bin(stan::io::var_context& context__,
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
        static const char* function__ = "model_wei_dex_bin_namespace::model_wei_dex_bin";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 23;
            context__.validate_dims("data initialization", "M", "int", context__.to_vec());
            M = int(0);
            vals_i__ = context__.vals_i("M");
            pos__ = 0;
            M = vals_i__[pos__++];
            check_greater_or_equal(function__, "M", M, 1);
            current_statement_begin__ = 24;
            context__.validate_dims("data initialization", "Mgen", "int", context__.to_vec());
            Mgen = int(0);
            vals_i__ = context__.vals_i("Mgen");
            pos__ = 0;
            Mgen = vals_i__[pos__++];
            check_greater_or_equal(function__, "Mgen", Mgen, 1);
            current_statement_begin__ = 25;
            context__.validate_dims("data initialization", "Nt", "int", context__.to_vec());
            Nt = int(0);
            vals_i__ = context__.vals_i("Nt");
            pos__ = 0;
            Nt = vals_i__[pos__++];
            check_greater_or_equal(function__, "Nt", Nt, 1);
            current_statement_begin__ = 26;
            validate_non_negative_index("y", "M", M);
            validate_non_negative_index("y", "Nt", Nt);
            context__.validate_dims("data initialization", "y", "double", context__.to_vec(M,Nt));
            y = std::vector<std::vector<double> >(M, std::vector<double>(Nt, double(0)));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_k_0_max__ = M;
            size_t y_k_1_max__ = Nt;
            for (size_t k_1__ = 0; k_1__ < y_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                    y[k_0__][k_1__] = vals_r__[pos__++];
                }
            }
            size_t y_i_0_max__ = M;
            size_t y_i_1_max__ = Nt;
            for (size_t i_0__ = 0; i_0__ < y_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < y_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "y[i_0__][i_1__]", y[i_0__][i_1__], 0);
                }
            }
            current_statement_begin__ = 27;
            validate_non_negative_index("N", "M", M);
            context__.validate_dims("data initialization", "N", "int", context__.to_vec(M));
            N = std::vector<int>(M, int(0));
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            size_t N_k_0_max__ = M;
            for (size_t k_0__ = 0; k_0__ < N_k_0_max__; ++k_0__) {
                N[k_0__] = vals_i__[pos__++];
            }
            size_t N_i_0_max__ = M;
            for (size_t i_0__ = 0; i_0__ < N_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "N[i_0__]", N[i_0__], 0);
            }
            current_statement_begin__ = 28;
            validate_non_negative_index("cexploc0prior", "2", 2);
            context__.validate_dims("data initialization", "cexploc0prior", "double", context__.to_vec(2));
            cexploc0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("cexploc0prior");
            pos__ = 0;
            size_t cexploc0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < cexploc0prior_k_0_max__; ++k_0__) {
                cexploc0prior[k_0__] = vals_r__[pos__++];
            }
            size_t cexploc0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < cexploc0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "cexploc0prior[i_0__]", cexploc0prior[i_0__], 0);
            }
            current_statement_begin__ = 29;
            validate_non_negative_index("wexploc0prior", "2", 2);
            context__.validate_dims("data initialization", "wexploc0prior", "double", context__.to_vec(2));
            wexploc0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("wexploc0prior");
            pos__ = 0;
            size_t wexploc0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < wexploc0prior_k_0_max__; ++k_0__) {
                wexploc0prior[k_0__] = vals_r__[pos__++];
            }
            size_t wexploc0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < wexploc0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "wexploc0prior[i_0__]", wexploc0prior[i_0__], 0);
            }
            current_statement_begin__ = 30;
            validate_non_negative_index("cexpsc0prior", "2", 2);
            context__.validate_dims("data initialization", "cexpsc0prior", "double", context__.to_vec(2));
            cexpsc0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("cexpsc0prior");
            pos__ = 0;
            size_t cexpsc0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < cexpsc0prior_k_0_max__; ++k_0__) {
                cexpsc0prior[k_0__] = vals_r__[pos__++];
            }
            size_t cexpsc0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < cexpsc0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "cexpsc0prior[i_0__]", cexpsc0prior[i_0__], 0);
            }
            current_statement_begin__ = 31;
            validate_non_negative_index("wexpsc0prior", "2", 2);
            context__.validate_dims("data initialization", "wexpsc0prior", "double", context__.to_vec(2));
            wexpsc0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("wexpsc0prior");
            pos__ = 0;
            size_t wexpsc0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < wexpsc0prior_k_0_max__; ++k_0__) {
                wexpsc0prior[k_0__] = vals_r__[pos__++];
            }
            size_t wexpsc0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < wexpsc0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "wexpsc0prior[i_0__]", wexpsc0prior[i_0__], 0);
            }
            current_statement_begin__ = 32;
            validate_non_negative_index("pn0prior", "2", 2);
            context__.validate_dims("data initialization", "pn0prior", "double", context__.to_vec(2));
            pn0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("pn0prior");
            pos__ = 0;
            size_t pn0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < pn0prior_k_0_max__; ++k_0__) {
                pn0prior[k_0__] = vals_r__[pos__++];
            }
            size_t pn0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < pn0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "pn0prior[i_0__]", pn0prior[i_0__], 0);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 35;
            num_params_r__ += 1;
            current_statement_begin__ = 36;
            num_params_r__ += 1;
            current_statement_begin__ = 37;
            validate_non_negative_index("C", "M", M);
            num_params_r__ += (1 * M);
            current_statement_begin__ = 38;
            validate_non_negative_index("w", "M", M);
            num_params_r__ += (1 * M);
            current_statement_begin__ = 39;
            num_params_r__ += 1;
            current_statement_begin__ = 40;
            num_params_r__ += 1;
            current_statement_begin__ = 41;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_wei_dex_bin() { }
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
        current_statement_begin__ = 35;
        if (!(context__.contains_r("cloc")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable cloc missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("cloc");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "cloc", "double", context__.to_vec());
        double cloc(0);
        cloc = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, cloc);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable cloc: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 36;
        if (!(context__.contains_r("wloc")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable wloc missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("wloc");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "wloc", "double", context__.to_vec());
        double wloc(0);
        wloc = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, wloc);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable wloc: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 37;
        if (!(context__.contains_r("C")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable C missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("C");
        pos__ = 0U;
        validate_non_negative_index("C", "M", M);
        context__.validate_dims("parameter initialization", "C", "double", context__.to_vec(M));
        std::vector<double> C(M, double(0));
        size_t C_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < C_k_0_max__; ++k_0__) {
            C[k_0__] = vals_r__[pos__++];
        }
        size_t C_i_0_max__ = M;
        for (size_t i_0__ = 0; i_0__ < C_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_lb_unconstrain(cloc, C[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable C: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 38;
        if (!(context__.contains_r("w")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable w missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("w");
        pos__ = 0U;
        validate_non_negative_index("w", "M", M);
        context__.validate_dims("parameter initialization", "w", "double", context__.to_vec(M));
        std::vector<double> w(M, double(0));
        size_t w_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < w_k_0_max__; ++k_0__) {
            w[k_0__] = vals_r__[pos__++];
        }
        size_t w_i_0_max__ = M;
        for (size_t i_0__ = 0; i_0__ < w_i_0_max__; ++i_0__) {
            try {
                writer__.scalar_lb_unconstrain(wloc, w[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable w: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        current_statement_begin__ = 39;
        if (!(context__.contains_r("csc")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable csc missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("csc");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "csc", "double", context__.to_vec());
        double csc(0);
        csc = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, csc);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable csc: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 40;
        if (!(context__.contains_r("wsc")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable wsc missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("wsc");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "wsc", "double", context__.to_vec());
        double wsc(0);
        wsc = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, wsc);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable wsc: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 41;
        if (!(context__.contains_r("pn")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable pn missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("pn");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "pn", "double", context__.to_vec());
        double pn(0);
        pn = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, pn);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable pn: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 35;
            local_scalar_t__ cloc;
            (void) cloc;  // dummy to suppress unused var warning
            if (jacobian__)
                cloc = in__.scalar_lb_constrain(0, lp__);
            else
                cloc = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 36;
            local_scalar_t__ wloc;
            (void) wloc;  // dummy to suppress unused var warning
            if (jacobian__)
                wloc = in__.scalar_lb_constrain(0, lp__);
            else
                wloc = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 37;
            std::vector<local_scalar_t__> C;
            size_t C_d_0_max__ = M;
            C.reserve(C_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < C_d_0_max__; ++d_0__) {
                if (jacobian__)
                    C.push_back(in__.scalar_lb_constrain(cloc, lp__));
                else
                    C.push_back(in__.scalar_lb_constrain(cloc));
            }
            current_statement_begin__ = 38;
            std::vector<local_scalar_t__> w;
            size_t w_d_0_max__ = M;
            w.reserve(w_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < w_d_0_max__; ++d_0__) {
                if (jacobian__)
                    w.push_back(in__.scalar_lb_constrain(wloc, lp__));
                else
                    w.push_back(in__.scalar_lb_constrain(wloc));
            }
            current_statement_begin__ = 39;
            local_scalar_t__ csc;
            (void) csc;  // dummy to suppress unused var warning
            if (jacobian__)
                csc = in__.scalar_lb_constrain(0, lp__);
            else
                csc = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 40;
            local_scalar_t__ wsc;
            (void) wsc;  // dummy to suppress unused var warning
            if (jacobian__)
                wsc = in__.scalar_lb_constrain(0, lp__);
            else
                wsc = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 41;
            local_scalar_t__ pn;
            (void) pn;  // dummy to suppress unused var warning
            if (jacobian__)
                pn = in__.scalar_lub_constrain(0, 1, lp__);
            else
                pn = in__.scalar_lub_constrain(0, 1);
            // model body
            current_statement_begin__ = 45;
            lp_accum__.add(beta_log<propto__>(pn, get_base1(pn0prior, 1, "pn0prior", 1), get_base1(pn0prior, 2, "pn0prior", 1)));
            current_statement_begin__ = 52;
            lp_accum__.add(gamma_log<propto__>(cloc, get_base1(cexploc0prior, 1, "cexploc0prior", 1), get_base1(cexploc0prior, 2, "cexploc0prior", 1)));
            current_statement_begin__ = 53;
            lp_accum__.add(gamma_log<propto__>(wloc, get_base1(wexploc0prior, 1, "wexploc0prior", 1), get_base1(wexploc0prior, 2, "wexploc0prior", 1)));
            current_statement_begin__ = 54;
            lp_accum__.add(gamma_log<propto__>(csc, get_base1(cexpsc0prior, 1, "cexpsc0prior", 1), get_base1(cexpsc0prior, 2, "cexpsc0prior", 1)));
            current_statement_begin__ = 55;
            lp_accum__.add(gamma_log<propto__>(wsc, get_base1(wexpsc0prior, 1, "wexpsc0prior", 1), get_base1(wexpsc0prior, 2, "wexpsc0prior", 1)));
            current_statement_begin__ = 58;
            for (int m = 1; m <= M; ++m) {
                current_statement_begin__ = 59;
                lp_accum__.add(exploc_lpdf<propto__>(get_base1(C, m, "C", 1), cloc, csc, pstream__));
                current_statement_begin__ = 60;
                lp_accum__.add(exploc_lpdf<propto__>(get_base1(w, m, "w", 1), wloc, wsc, pstream__));
                current_statement_begin__ = 61;
                lp_accum__.add(binomial_log(get_base1(N, m, "N", 1), Nt, pn));
                current_statement_begin__ = 62;
                for (int j = 1; j <= Nt; ++j) {
                    current_statement_begin__ = 63;
                    if (as_bool(logical_gt(get_base1(get_base1(y, m, "y", 1), j, "y", 2), 1e-6))) {
                        current_statement_begin__ = 64;
                        lp_accum__.add(weibull_log<propto__>(get_base1(get_base1(y, m, "y", 1), j, "y", 2), get_base1(w, m, "w", 1), get_base1(C, m, "C", 1)));
                    }
                }
            }
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
        names__.push_back("cloc");
        names__.push_back("wloc");
        names__.push_back("C");
        names__.push_back("w");
        names__.push_back("csc");
        names__.push_back("wsc");
        names__.push_back("pn");
        names__.push_back("Ngen");
        names__.push_back("Cgen");
        names__.push_back("wgen");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(M);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(M);
        dimss__.push_back(dims__);
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
        dims__.push_back(Mgen);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(Mgen);
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
        static const char* function__ = "model_wei_dex_bin_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double cloc = in__.scalar_lb_constrain(0);
        vars__.push_back(cloc);
        double wloc = in__.scalar_lb_constrain(0);
        vars__.push_back(wloc);
        std::vector<double> C;
        size_t C_d_0_max__ = M;
        C.reserve(C_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < C_d_0_max__; ++d_0__) {
            C.push_back(in__.scalar_lb_constrain(cloc));
        }
        size_t C_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < C_k_0_max__; ++k_0__) {
            vars__.push_back(C[k_0__]);
        }
        std::vector<double> w;
        size_t w_d_0_max__ = M;
        w.reserve(w_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < w_d_0_max__; ++d_0__) {
            w.push_back(in__.scalar_lb_constrain(wloc));
        }
        size_t w_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < w_k_0_max__; ++k_0__) {
            vars__.push_back(w[k_0__]);
        }
        double csc = in__.scalar_lb_constrain(0);
        vars__.push_back(csc);
        double wsc = in__.scalar_lb_constrain(0);
        vars__.push_back(wsc);
        double pn = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(pn);
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
            validate_non_negative_index("Ngen", "Mgen", Mgen);
            std::vector<int> Ngen(Mgen, int(0));
            stan::math::fill(Ngen, std::numeric_limits<int>::min());
            current_statement_begin__ = 75;
            validate_non_negative_index("Cgen", "Mgen", Mgen);
            std::vector<double> Cgen(Mgen, double(0));
            stan::math::initialize(Cgen, DUMMY_VAR__);
            stan::math::fill(Cgen, DUMMY_VAR__);
            current_statement_begin__ = 76;
            validate_non_negative_index("wgen", "Mgen", Mgen);
            std::vector<double> wgen(Mgen, double(0));
            stan::math::initialize(wgen, DUMMY_VAR__);
            stan::math::fill(wgen, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 90;
            for (int m = 1; m <= Mgen; ++m) {
                current_statement_begin__ = 95;
                stan::model::assign(Cgen, 
                            stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                            (cloc + exponential_rng(inv(csc), base_rng__)), 
                            "assigning variable Cgen");
                current_statement_begin__ = 96;
                stan::model::assign(wgen, 
                            stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                            (wloc + exponential_rng(inv(wsc), base_rng__)), 
                            "assigning variable wgen");
                current_statement_begin__ = 97;
                stan::model::assign(Ngen, 
                            stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                            binomial_rng(Nt, pn, base_rng__), 
                            "assigning variable Ngen");
            }
            // validate, write generated quantities
            current_statement_begin__ = 74;
            size_t Ngen_i_0_max__ = Mgen;
            for (size_t i_0__ = 0; i_0__ < Ngen_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "Ngen[i_0__]", Ngen[i_0__], 0);
            }
            size_t Ngen_k_0_max__ = Mgen;
            for (size_t k_0__ = 0; k_0__ < Ngen_k_0_max__; ++k_0__) {
                vars__.push_back(Ngen[k_0__]);
            }
            current_statement_begin__ = 75;
            size_t Cgen_i_0_max__ = Mgen;
            for (size_t i_0__ = 0; i_0__ < Cgen_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "Cgen[i_0__]", Cgen[i_0__], 0);
            }
            size_t Cgen_k_0_max__ = Mgen;
            for (size_t k_0__ = 0; k_0__ < Cgen_k_0_max__; ++k_0__) {
                vars__.push_back(Cgen[k_0__]);
            }
            current_statement_begin__ = 76;
            size_t wgen_i_0_max__ = Mgen;
            for (size_t i_0__ = 0; i_0__ < wgen_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "wgen[i_0__]", wgen[i_0__], 0);
            }
            size_t wgen_k_0_max__ = Mgen;
            for (size_t k_0__ = 0; k_0__ < wgen_k_0_max__; ++k_0__) {
                vars__.push_back(wgen[k_0__]);
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
        return "model_wei_dex_bin";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "cloc";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "wloc";
        param_names__.push_back(param_name_stream__.str());
        size_t C_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < C_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "C" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t w_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < w_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "w" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "csc";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "wsc";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pn";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t Ngen_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < Ngen_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Ngen" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t Cgen_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < Cgen_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Cgen" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t wgen_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < wgen_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "wgen" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "cloc";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "wloc";
        param_names__.push_back(param_name_stream__.str());
        size_t C_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < C_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "C" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t w_k_0_max__ = M;
        for (size_t k_0__ = 0; k_0__ < w_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "w" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "csc";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "wsc";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pn";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t Ngen_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < Ngen_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Ngen" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t Cgen_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < Cgen_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Cgen" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t wgen_k_0_max__ = Mgen;
        for (size_t k_0__ = 0; k_0__ < wgen_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "wgen" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_wei_dex_bin_namespace::model_wei_dex_bin stan_model;
#endif
