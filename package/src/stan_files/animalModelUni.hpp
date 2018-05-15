/*
    stanAnimal is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    stanAnimal is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with stanAnimal.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.17.0

#include <stan/model/model_header.hpp>

namespace model_animalModelUni_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_animalModelUni");
    reader.add_event(42, 42, "end", "model_animalModelUni");
    return reader;
}

#include <meta_header.hpp>
 class model_animalModelUni : public prob_grad {
private:
    int J;
    int N;
    vector<vector_d> X;
    vector<double> Y;
    matrix_d A;
    matrix_d LA;
public:
    model_animalModelUni(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_animalModelUni(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_animalModelUni_namespace::model_animalModelUni";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "J", "int", context__.to_vec());
            J = int(0);
            vals_i__ = context__.vals_i("J");
            pos__ = 0;
            J = vals_i__[pos__++];
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 4;
            validate_non_negative_index("X", "N", N);
            validate_non_negative_index("X", "J", J);
            context__.validate_dims("data initialization", "X", "vector_d", context__.to_vec(N,J));
            validate_non_negative_index("X", "N", N);
            validate_non_negative_index("X", "J", J);
            X = std::vector<vector_d>(N,vector_d(static_cast<Eigen::VectorXd::Index>(J)));
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_i_vec_lim__ = J;
            for (size_t i_vec__ = 0; i_vec__ < X_i_vec_lim__; ++i_vec__) {
                size_t X_limit_0__ = N;
                for (size_t i_0__ = 0; i_0__ < X_limit_0__; ++i_0__) {
                    X[i_0__][i_vec__] = vals_r__[pos__++];
            }
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("Y", "N", N);
            context__.validate_dims("data initialization", "Y", "double", context__.to_vec(N));
            validate_non_negative_index("Y", "N", N);
            Y = std::vector<double>(N,double(0));
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < Y_limit_0__; ++i_0__) {
                Y[i_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("A", "N", N);
            validate_non_negative_index("A", "N", N);
            context__.validate_dims("data initialization", "A", "matrix_d", context__.to_vec(N,N));
            validate_non_negative_index("A", "N", N);
            validate_non_negative_index("A", "N", N);
            A = matrix_d(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
            vals_r__ = context__.vals_r("A");
            pos__ = 0;
            size_t A_k_mat_lim__ = N;
            for (size_t n_mat__ = 0; n_mat__ < A_k_mat_lim__; ++n_mat__) {
                for (size_t m_mat__ = 0; m_mat__ < A_k_mat_lim__; ++m_mat__) {
                    A(m_mat__,n_mat__) = vals_r__[pos__++];
                }
            }

            // validate, data variables
            current_statement_begin__ = 2;
            check_greater_or_equal(function__,"J",J,1);
            current_statement_begin__ = 3;
            check_greater_or_equal(function__,"N",N,0);
            current_statement_begin__ = 4;
            current_statement_begin__ = 5;
            current_statement_begin__ = 6;
            stan::math::check_cov_matrix(function__,"A",A);
            // initialize data variables
            current_statement_begin__ = 9;
            validate_non_negative_index("LA", "N", N);
            validate_non_negative_index("LA", "N", N);
            LA = matrix_d(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
            stan::math::fill(LA,DUMMY_VAR__);

            current_statement_begin__ = 10;
            stan::math::assign(LA, cholesky_decompose(A));

            // validate transformed data
            current_statement_begin__ = 9;

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 13;
            validate_non_negative_index("a_tilde", "N", N);
            num_params_r__ += N;
            current_statement_begin__ = 14;
        validate_non_negative_index("beta", "J", J);
            num_params_r__ += J;
            current_statement_begin__ = 17;
            ++num_params_r__;
            current_statement_begin__ = 20;
            ++num_params_r__;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_animalModelUni() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("a_tilde")))
            throw std::runtime_error("variable a_tilde missing");
        vals_r__ = context__.vals_r("a_tilde");
        pos__ = 0U;
        validate_non_negative_index("a_tilde", "N", N);
        context__.validate_dims("initialization", "a_tilde", "vector_d", context__.to_vec(N));
        vector_d a_tilde(static_cast<Eigen::VectorXd::Index>(N));
        for (int j1__ = 0U; j1__ < N; ++j1__)
            a_tilde(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_unconstrain(a_tilde);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable a_tilde: ") + e.what());
        }

        if (!(context__.contains_r("beta")))
            throw std::runtime_error("variable beta missing");
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "J", J);
        context__.validate_dims("initialization", "beta", "row_vector_d", context__.to_vec(J));
        row_vector_d beta(static_cast<Eigen::VectorXd::Index>(J));
        for (int j1__ = 0U; j1__ < J; ++j1__)
            beta(j1__) = vals_r__[pos__++];
        try {
            writer__.row_vector_unconstrain(beta);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable beta: ") + e.what());
        }

        if (!(context__.contains_r("sigma_G")))
            throw std::runtime_error("variable sigma_G missing");
        vals_r__ = context__.vals_r("sigma_G");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigma_G", "double", context__.to_vec());
        double sigma_G(0);
        sigma_G = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,sigma_G);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma_G: ") + e.what());
        }

        if (!(context__.contains_r("sigma_R")))
            throw std::runtime_error("variable sigma_R missing");
        vals_r__ = context__.vals_r("sigma_R");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigma_R", "double", context__.to_vec());
        double sigma_R(0);
        sigma_R = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,sigma_R);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma_R: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<T__> in__(params_r__,params_i__);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  a_tilde;
            (void) a_tilde;  // dummy to suppress unused var warning
            if (jacobian__)
                a_tilde = in__.vector_constrain(N,lp__);
            else
                a_tilde = in__.vector_constrain(N);

            Eigen::Matrix<T__,1,Eigen::Dynamic>  beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.row_vector_constrain(J,lp__);
            else
                beta = in__.row_vector_constrain(J);

            T__ sigma_G;
            (void) sigma_G;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_G = in__.scalar_lb_constrain(0,lp__);
            else
                sigma_G = in__.scalar_lb_constrain(0);

            T__ sigma_R;
            (void) sigma_R;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_R = in__.scalar_lb_constrain(0,lp__);
            else
                sigma_R = in__.scalar_lb_constrain(0);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body
            {
            current_statement_begin__ = 23;
            validate_non_negative_index("mu", "N", N);
            Eigen::Matrix<T__,Eigen::Dynamic,1>  mu(static_cast<Eigen::VectorXd::Index>(N));
            (void) mu;  // dummy to suppress unused var warning

            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu,DUMMY_VAR__);
            current_statement_begin__ = 24;
            validate_non_negative_index("a", "N", N);
            Eigen::Matrix<T__,Eigen::Dynamic,1>  a(static_cast<Eigen::VectorXd::Index>(N));
            (void) a;  // dummy to suppress unused var warning

            stan::math::initialize(a, DUMMY_VAR__);
            stan::math::fill(a,DUMMY_VAR__);


            current_statement_begin__ = 26;
            lp_accum__.add(normal_log<propto__>(a_tilde, 0, 1));
            current_statement_begin__ = 27;
            stan::math::assign(a, multiply(sqrt(sigma_G),multiply(LA,a_tilde)));
            current_statement_begin__ = 29;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 30;
                stan::math::assign(get_base1_lhs(mu,n,"mu",1), (multiply(beta,get_base1(X,n,"X",1)) + get_base1(a,n,"a",1)));
            }
            current_statement_begin__ = 32;
            lp_accum__.add(normal_log<propto__>(Y, mu, sigma_R));
            current_statement_begin__ = 34;
            lp_accum__.add(normal_log<propto__>(to_vector(beta), 0, 1));
            current_statement_begin__ = 36;
            lp_accum__.add(normal_log<propto__>(sigma_G, 0, 1));
            current_statement_begin__ = 37;
            lp_accum__.add(normal_log<propto__>(sigma_R, 0, 1));
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
        names__.push_back("a_tilde");
        names__.push_back("beta");
        names__.push_back("sigma_G");
        names__.push_back("sigma_R");
        names__.push_back("sigma_E");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(J);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
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
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "model_animalModelUni_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d a_tilde = in__.vector_constrain(N);
        row_vector_d beta = in__.row_vector_constrain(J);
        double sigma_G = in__.scalar_lb_constrain(0);
        double sigma_R = in__.scalar_lb_constrain(0);
            for (int k_0__ = 0; k_0__ < N; ++k_0__) {
            vars__.push_back(a_tilde[k_0__]);
            }
            for (int k_0__ = 0; k_0__ < J; ++k_0__) {
            vars__.push_back(beta[k_0__]);
            }
        vars__.push_back(sigma_G);
        vars__.push_back(sigma_R);

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters

            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 40;
            double sigma_E(0.0);
            (void) sigma_E;  // dummy to suppress unused var warning

            stan::math::initialize(sigma_E, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(sigma_E,DUMMY_VAR__);


            current_statement_begin__ = 41;
            stan::math::assign(sigma_E, (sigma_R * sigma_R));

            // validate generated quantities
            current_statement_begin__ = 40;

            // write generated quantities
        vars__.push_back(sigma_E);

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
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_animalModelUni";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "a_tilde" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= J; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_G";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_R";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_E";
        param_names__.push_back(param_name_stream__.str());
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "a_tilde" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= J; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_G";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_R";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_E";
        param_names__.push_back(param_name_stream__.str());
    }

}; // model

}

typedef model_animalModelUni_namespace::model_animalModelUni stan_model;


#endif