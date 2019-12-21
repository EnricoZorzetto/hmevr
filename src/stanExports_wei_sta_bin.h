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
namespace model_wei_sta_bin_namespace {
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
    reader.add_event(0, 0, "start", "model_wei_sta_bin");
    reader.add_event(88, 86, "end", "model_wei_sta_bin");
    return reader;
}
#include <stan_meta_header.hpp>
class model_wei_sta_bin : public prob_grad {
private:
        int M;
        int Mgen;
        int Nt;
        std::vector<std::vector<double> > y;
        std::vector<int> N;
        std::vector<double> C0prior;
        std::vector<double> w0prior;
        std::vector<double> pn0prior;
public:
    model_wei_sta_bin(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_wei_sta_bin(stan::io::var_context& context__,
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
        static const char* function__ = "model_wei_sta_bin_namespace::model_wei_sta_bin";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "M", "int", context__.to_vec());
            M = int(0);
            vals_i__ = context__.vals_i("M");
            pos__ = 0;
            M = vals_i__[pos__++];
            check_greater_or_equal(function__, "M", M, 1);
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "Mgen", "int", context__.to_vec());
            Mgen = int(0);
            vals_i__ = context__.vals_i("Mgen");
            pos__ = 0;
            Mgen = vals_i__[pos__++];
            check_greater_or_equal(function__, "Mgen", Mgen, 1);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "Nt", "int", context__.to_vec());
            Nt = int(0);
            vals_i__ = context__.vals_i("Nt");
            pos__ = 0;
            Nt = vals_i__[pos__++];
            check_greater_or_equal(function__, "Nt", Nt, 1);
            current_statement_begin__ = 5;
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
            current_statement_begin__ = 6;
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
            current_statement_begin__ = 9;
            validate_non_negative_index("C0prior", "2", 2);
            context__.validate_dims("data initialization", "C0prior", "double", context__.to_vec(2));
            C0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("C0prior");
            pos__ = 0;
            size_t C0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < C0prior_k_0_max__; ++k_0__) {
                C0prior[k_0__] = vals_r__[pos__++];
            }
            size_t C0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < C0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "C0prior[i_0__]", C0prior[i_0__], 0);
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("w0prior", "2", 2);
            context__.validate_dims("data initialization", "w0prior", "double", context__.to_vec(2));
            w0prior = std::vector<double>(2, double(0));
            vals_r__ = context__.vals_r("w0prior");
            pos__ = 0;
            size_t w0prior_k_0_max__ = 2;
            for (size_t k_0__ = 0; k_0__ < w0prior_k_0_max__; ++k_0__) {
                w0prior[k_0__] = vals_r__[pos__++];
            }
            size_t w0prior_i_0_max__ = 2;
            for (size_t i_0__ = 0; i_0__ < w0prior_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "w0prior[i_0__]", w0prior[i_0__], 0);
            }
            current_statement_begin__ = 11;
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
            current_statement_begin__ = 18;
            num_params_r__ += 1;
            current_statement_begin__ = 19;
            num_params_r__ += 1;
            current_statement_begin__ = 20;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_wei_sta_bin() { }
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
        current_statement_begin__ = 18;
        if (!(context__.contains_r("C")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable C missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("C");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "C", "double", context__.to_vec());
        double C(0);
        C = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, C);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable C: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 19;
        if (!(context__.contains_r("w")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable w missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("w");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "w", "double", context__.to_vec());
        double w(0);
        w = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, w);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable w: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 20;
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
            current_statement_begin__ = 18;
            local_scalar_t__ C;
            (void) C;  // dummy to suppress unused var warning
            if (jacobian__)
                C = in__.scalar_lb_constrain(0, lp__);
            else
                C = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 19;
            local_scalar_t__ w;
            (void) w;  // dummy to suppress unused var warning
            if (jacobian__)
                w = in__.scalar_lb_constrain(0, lp__);
            else
                w = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 20;
            local_scalar_t__ pn;
            (void) pn;  // dummy to suppress unused var warning
            if (jacobian__)
                pn = in__.scalar_lub_constrain(0, 1, lp__);
            else
                pn = in__.scalar_lub_constrain(0, 1);
            // model body
            current_statement_begin__ = 24;
            lp_accum__.add(beta_log<propto__>(pn, get_base1(pn0prior, 1, "pn0prior", 1), get_base1(pn0prior, 2, "pn0prior", 1)));
            current_statement_begin__ = 25;
            lp_accum__.add(gamma_log<propto__>(C, get_base1(C0prior, 1, "C0prior", 1), get_base1(C0prior, 2, "C0prior", 1)));
            current_statement_begin__ = 26;
            lp_accum__.add(gamma_log<propto__>(w, get_base1(w0prior, 1, "w0prior", 1), get_base1(w0prior, 2, "w0prior", 1)));
            current_statement_begin__ = 28;
            for (int m = 1; m <= M; ++m) {
                current_statement_begin__ = 29;
                lp_accum__.add(binomial_log(get_base1(N, m, "N", 1), Nt, pn));
                current_statement_begin__ = 30;
                for (int j = 1; j <= Nt; ++j) {
                    current_statement_begin__ = 31;
                    if (as_bool(logical_gt(get_base1(get_base1(y, m, "y", 1), j, "y", 2), 1e-6))) {
                        current_statement_begin__ = 32;
                        lp_accum__.add(weibull_log(get_base1(get_base1(y, m, "y", 1), j, "y", 2), w, C));
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
        names__.push_back("C");
        names__.push_back("w");
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
        static const char* function__ = "model_wei_sta_bin_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double C = in__.scalar_lb_constrain(0);
        vars__.push_back(C);
        double w = in__.scalar_lb_constrain(0);
        vars__.push_back(w);
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
            current_statement_begin__ = 42;
            validate_non_negative_index("Ngen", "Mgen", Mgen);
            std::vector<int> Ngen(Mgen, int(0));
            stan::math::fill(Ngen, std::numeric_limits<int>::min());
            current_statement_begin__ = 43;
            validate_non_negative_index("Cgen", "Mgen", Mgen);
            std::vector<double> Cgen(Mgen, double(0));
            stan::math::initialize(Cgen, DUMMY_VAR__);
            stan::math::fill(Cgen, DUMMY_VAR__);
            current_statement_begin__ = 44;
            validate_non_negative_index("wgen", "Mgen", Mgen);
            std::vector<double> wgen(Mgen, double(0));
            stan::math::initialize(wgen, DUMMY_VAR__);
            stan::math::fill(wgen, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 58;
            for (int m = 1; m <= Mgen; ++m) {
                current_statement_begin__ = 59;
                stan::model::assign(Cgen, 
                            stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                            C, 
                            "assigning variable Cgen");
                current_statement_begin__ = 60;
                stan::model::assign(wgen, 
                            stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                            w, 
                            "assigning variable wgen");
                current_statement_begin__ = 61;
                stan::model::assign(Ngen, 
                            stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                            binomial_rng(Nt, pn, base_rng__), 
                            "assigning variable Ngen");
            }
            // validate, write generated quantities
            current_statement_begin__ = 42;
            size_t Ngen_i_0_max__ = Mgen;
            for (size_t i_0__ = 0; i_0__ < Ngen_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "Ngen[i_0__]", Ngen[i_0__], 0);
            }
            size_t Ngen_k_0_max__ = Mgen;
            for (size_t k_0__ = 0; k_0__ < Ngen_k_0_max__; ++k_0__) {
                vars__.push_back(Ngen[k_0__]);
            }
            current_statement_begin__ = 43;
            size_t Cgen_i_0_max__ = Mgen;
            for (size_t i_0__ = 0; i_0__ < Cgen_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "Cgen[i_0__]", Cgen[i_0__], 0);
            }
            size_t Cgen_k_0_max__ = Mgen;
            for (size_t k_0__ = 0; k_0__ < Cgen_k_0_max__; ++k_0__) {
                vars__.push_back(Cgen[k_0__]);
            }
            current_statement_begin__ = 44;
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
        return "model_wei_sta_bin";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "C";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "w";
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
        param_name_stream__ << "C";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "w";
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
typedef model_wei_sta_bin_namespace::model_wei_sta_bin stan_model;
#endif