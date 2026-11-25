/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_OPTIMIZER_CMAES_HH
#define OPTIMIST_OPTIMIZER_CMAES_HH

#include <random>

#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the Covariance Matrix Adaptation Evolution
     * Strategy (CMA-ES), with optional IPOP restarts.
     *
     * \includedoc docs/markdown/Optimizer/CMAES.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class CMAES : public Optimizer<Vector, CMAES<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using DenseMatrix = Utils::DenseMatrix<Vector>;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_sigma{0};        /**< Initial step size (0 means automatic). */
      Integer m_population{0};  /**< Population size (0 means automatic). */
      Integer m_restarts{0};    /**< Number of IPOP restarts. */
      std::mt19937 m_generator{42}; /**< Random number generator. */

     public:
      /**
       * Class constructor for the CMA-ES solver.
       */
      CMAES() {}

      /**
       * Get the CMA-ES solver name.
       */
      constexpr std::string name_impl() const {
        return "CMA-ES";
      }

      /**
       * Set the seed of the random number generator.
       */
      void seed(const unsigned t_seed) {
        this->m_generator.seed(t_seed);
      }

      /**
       * Get the initial step size (0 means automatic: 0.3 times the bounded
       * box width, or 1 if the bounds are not finite).
       */
      Scalar sigma() const {
        return this->m_sigma;
      }

      /**
       * Set the initial step size.
       */
      void sigma(const Scalar t_sigma) {
        OPTIMIST_ASSERT(t_sigma >= 0,
                        "Optimist::Optimizer::CMAES::sigma(...): invalid "
                        "input detected.");
        this->m_sigma = t_sigma;
      }

      /**
       * Get the population size (0 means \f$ 4 + \lfloor 3\ln n \rfloor \f$).
       */
      Integer population() const {
        return this->m_population;
      }

      /**
       * Set the population size.
       */
      void population(const Integer t_population) {
        OPTIMIST_ASSERT(t_population == 0 || t_population >= 4,
                        "Optimist::Optimizer::CMAES::population(...): invalid "
                        "input detected.");
        this->m_population = t_population;
      }

      /**
       * Get the number of IPOP restarts.
       */
      Integer restarts() const {
        return this->m_restarts;
      }

      /**
       * Set the number of IPOP restarts: each restart doubles the population
       * size and starts from a random point of the bounded box (or from the
       * initial guess if the bounds are not finite).
       */
      void restarts(const Integer t_restarts) {
        OPTIMIST_ASSERT(t_restarts >= 0,
                        "Optimist::Optimizer::CMAES::restarts(...): invalid "
                        "input detected.");
        this->m_restarts = t_restarts;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] x_ini Initialization point (initial mean).
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda &&function,
                      const Vector &x_ini,
                      Vector &x_sol) {
        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        const Integer n{static_cast<Integer>(x_ini.size())};
        const bool bounded{this->has_finite_bounds()};
        Dense lb(n), ub(n);
        if (bounded) {
          lb = Utils::to_dense(this->m_lower_bound);
          ub = Utils::to_dense(this->m_upper_bound);
        }
        const Scalar sigma_0{this->m_sigma > 0 ? this->m_sigma
                             : bounded ? Scalar(0.3) * (ub - lb).maxCoeff()
                                       : Scalar(1)};
        Integer lambda{this->m_population > 0
                           ? this->m_population
                           : 4 + static_cast<Integer>(
                                     3 * std::log(static_cast<Scalar>(n)))};

        Dense x_best(Utils::to_dense(x_ini));
        Scalar f_best{INFTY};
        bool converged_best{false};
        Dense mean(x_best);
        std::uniform_real_distribution<Scalar> unif(0, 1);

        this->m_iterations = 0;
        for (Integer run{0}; run <= this->m_restarts; ++run) {
          if (run > 0) {
            lambda *= 2;
            if (bounded) {
              for (Integer j{0}; j < n; ++j) {
                mean(j) = lb(j) + unif(this->m_generator) * (ub(j) - lb(j));
              }
            } else {
              mean = Utils::to_dense(x_ini);
            }
          }
          Dense x_run(n);
          Scalar f_run{INFTY};
          const bool converged{this->run(function,
                                         mean,
                                         sigma_0,
                                         lambda,
                                         bounded,
                                         lb,
                                         ub,
                                         x_run,
                                         f_run)};
          if (f_run < f_best) {
            f_best         = f_run;
            x_best         = x_run;
            converged_best = converged;
          }
          if (this->m_iterations >= this->m_max_iterations ||
              this->m_function_evaluations + lambda >
                  this->m_max_function_evaluations) {
            break;
          }
        }

        this->m_converged = converged_best;
        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = Utils::from_dense<Vector>(x_best);
        return this->m_converged;
      }

     private:
      template <typename FunctionLambda>
      Scalar eval(FunctionLambda &function, const Dense &x) {
        Scalar f;
        if (!this->evaluate_function(function, Utils::from_dense<Vector>(x), f) ||
            !std::isfinite(f)) {
          return INFTY;
        }
        return f;
      }

      /**
       * Single CMA-ES run (Hansen's tutorial, 2016).
       * \return True if the run met the tolerance on the step size.
       */
      template <typename FunctionLambda>
      bool run(FunctionLambda &function,
               Dense mean,
               Scalar sigma,
               const Integer lambda,
               const bool bounded,
               const Dense &lb,
               const Dense &ub,
               Dense &x_best,
               Scalar &f_best) {
        const Integer n{static_cast<Integer>(mean.size())};
        const Scalar sn{static_cast<Scalar>(n)};
        const Integer mu{lambda / 2};

        // Recombination weights
        Dense w(mu);
        for (Integer i{0}; i < mu; ++i) {
          w(i) = std::log(static_cast<Scalar>(mu) + Scalar(0.5)) -
                 std::log(static_cast<Scalar>(i + 1));
        }
        w /= w.sum();
        const Scalar mu_eff{1 / w.squaredNorm()};

        // Adaptation constants
        const Scalar c_s{(mu_eff + 2) / (sn + mu_eff + 5)};
        const Scalar d_s{
          1 + 2 * std::max(Scalar(0), std::sqrt((mu_eff - 1) / (sn + 1)) - 1) +
          c_s};
        const Scalar c_c{(4 + mu_eff / sn) / (sn + 4 + 2 * mu_eff / sn)};
        const Scalar c_1{2 / ((sn + Scalar(1.3)) * (sn + Scalar(1.3)) + mu_eff)};
        const Scalar c_mu{std::min(
            1 - c_1,
            2 * (mu_eff - 2 + 1 / mu_eff) / ((sn + 2) * (sn + 2) + mu_eff))};
        const Scalar chi_n{std::sqrt(sn) *
                           (1 - 1 / (4 * sn) + 1 / (21 * sn * sn))};

        // State
        Dense p_s(Dense::Zero(n)), p_c(Dense::Zero(n)), D(Dense::Ones(n));
        DenseMatrix C(DenseMatrix::Identity(n, n)),
            B(DenseMatrix::Identity(n, n));
        DenseMatrix X(n, lambda), Y(n, lambda);
        Dense f(lambda);
        std::vector<Integer> order(lambda);
        std::normal_distribution<Scalar> normal(0, 1);
        f_best = INFTY;
        x_best = mean;

        for (Integer generation{1};
             this->m_iterations < this->m_max_iterations;
             ++generation) {
          ++this->m_iterations;
          if (this->m_function_evaluations + lambda >
              this->m_max_function_evaluations) {
            return false;
          }

          // Sample and evaluate the population
          for (Integer k{0}; k < lambda; ++k) {
            Dense z(n);
            for (Integer j{0}; j < n; ++j) {
              z(j) = normal(this->m_generator);
            }
            X.col(k) = mean + sigma * (B * D.cwiseProduct(z));
            if (bounded) {
              X.col(k) = X.col(k).cwiseMax(lb).cwiseMin(ub);
            }
            Y.col(k) = (X.col(k) - mean) / sigma;
            f(k)     = this->eval(function, X.col(k));
          }
          std::iota(order.begin(), order.end(), Integer(0));
          std::sort(order.begin(),
                    order.end(),
                    [&f](const Integer a, const Integer b) {
                      return f(a) < f(b);
                    });
          if (f(order[0]) < f_best) {
            f_best = f(order[0]);
            x_best = X.col(order[0]);
          }

          // Mean update
          Dense y_w(Dense::Zero(n));
          for (Integer i{0}; i < mu; ++i) {
            y_w += w(i) * Y.col(order[i]);
          }
          mean += sigma * y_w;

          // Evolution paths
          const Dense C_inv_sqrt_y{B * (B.transpose() * y_w).cwiseQuotient(D)};
          p_s = (1 - c_s) * p_s + std::sqrt(c_s * (2 - c_s) * mu_eff) * C_inv_sqrt_y;
          const Scalar p_s_norm{p_s.norm()};
          const bool h_s{
            p_s_norm / std::sqrt(1 - std::pow(1 - c_s, Scalar(2 * generation))) <
            (Scalar(1.4) + 2 / (sn + 1)) * chi_n};
          p_c = (1 - c_c) * p_c +
                (h_s ? std::sqrt(c_c * (2 - c_c) * mu_eff) : Scalar(0)) * y_w;

          // Covariance matrix update
          DenseMatrix rank_mu(DenseMatrix::Zero(n, n));
          for (Integer i{0}; i < mu; ++i) {
            rank_mu += w(i) * Y.col(order[i]) * Y.col(order[i]).transpose();
          }
          const Scalar delta_h{h_s ? Scalar(0) : c_c * (2 - c_c)};
          C = (1 - c_1 - c_mu) * C +
              c_1 * (p_c * p_c.transpose() + delta_h * C) + c_mu * rank_mu;
          C = (C + C.transpose()).eval() / 2;

          // Step size update
          sigma *= std::exp((c_s / d_s) * (p_s_norm / chi_n - 1));

          // Eigendecomposition C = B diag(D^2) B'
          const Eigen::SelfAdjointEigenSolver<DenseMatrix> eig(C);
          if (eig.info() != Eigen::Success) {
            return false;
          }
          B = eig.eigenvectors();
          D = eig.eigenvalues().cwiseMax(EPSILON).cwiseSqrt();

          // Stopping criteria
          const Scalar step{sigma * D.maxCoeff()};
          if (this->m_verbose) {
            this->info(step, "Run with lambda = " + std::to_string(lambda));
          }
          // TolX: the search distribution collapsed; TolFun: the sampled
          // values are all equal to rounding level, so the selection carries
          // no more information (with a step already small)
          const Scalar scale{1 + mean.norm()};
          const Scalar f_range{f(order.back()) - f(order.front())};
          const bool flat{f_range <= 100 * EPSILON *
                                         std::max(Scalar(1),
                                                  std::abs(f(order.front())))};
          if (step < this->m_tolerance * scale ||
              (flat && step < std::sqrt(this->m_tolerance) * scale)) {
            // Report the mean if it is better than the best sample
            if (this->m_function_evaluations <
                this->m_max_function_evaluations) {
              const Scalar f_mean{this->eval(function, mean)};
              if (f_mean < f_best) {
                f_best = f_mean;
                x_best = mean;
              }
            }
            return true;
          }
          if (!std::isfinite(step) || D.maxCoeff() > 1.0e7 * D.minCoeff()) {
            return false;
          }
        }
        return false;
      }

    };  // class CMAES

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_CMAES_HH
