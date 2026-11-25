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

#ifndef OPTIMIST_OPTIMIZER_ADAM_HH
#define OPTIMIST_OPTIMIZER_ADAM_HH

#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the Adam method (adaptive moment
     * estimation), with the optional AMSGrad correction and step decay.
     *
     * \includedoc docs/markdown/Optimizer/Adam.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class Adam : public Optimizer<Vector, Adam<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;
      using typename Optimizer<Vector, Adam<Vector>>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_learning_rate{1.0e-3}; /**< Step size \f$ \eta \f$. */
      Scalar m_beta_1{0.9};   /**< First moment decay \f$ \beta_1 \f$. */
      Scalar m_beta_2{0.999}; /**< Second moment decay \f$ \beta_2 \f$. */
      Scalar m_epsilon{1.0e-8}; /**< Denominator regularization. */
      Scalar m_decay{0.0}; /**< Step decay: \f$ \eta_k = \eta/(1 + \gamma k)
                              \f$. */
      bool m_amsgrad{false}; /**< AMSGrad flag. */

     public:
      /**
       * Class constructor for the Adam solver.
       */
      Adam() {
        this->m_max_iterations = 10000;
        this->m_max_function_evaluations       = 10001;
        this->m_max_first_derivative_evaluations = 10001;
      }

      /**
       * Get the Adam solver name.
       */
      constexpr std::string name_impl() const {
        return this->m_amsgrad ? "AMSGrad" : "Adam";
      }

      /**
       * Get the learning rate.
       */
      Scalar learning_rate() const {
        return this->m_learning_rate;
      }

      /**
       * Set the learning rate \f$ \eta > 0 \f$.
       */
      void learning_rate(const Scalar t_learning_rate) {
        OPTIMIST_ASSERT(t_learning_rate > 0,
                        "Optimist::Optimizer::Adam::learning_rate(...): "
                        "invalid input detected.");
        this->m_learning_rate = t_learning_rate;
      }

      /**
       * Set the moment decay rates \f$ 0 \leq \beta_1, \beta_2 < 1 \f$.
       */
      void betas(const Scalar t_beta_1, const Scalar t_beta_2) {
        OPTIMIST_ASSERT(t_beta_1 >= 0 && t_beta_1 < 1 && t_beta_2 >= 0 &&
                            t_beta_2 < 1,
                        "Optimist::Optimizer::Adam::betas(...): invalid input "
                        "detected.");
        this->m_beta_1 = t_beta_1;
        this->m_beta_2 = t_beta_2;
      }

      /**
       * Set the step decay \f$ \gamma \geq 0 \f$, so that \f$ \eta_k =
       * \eta/(1 + \gamma k) \f$.
       */
      void decay(const Scalar t_decay) {
        OPTIMIST_ASSERT(t_decay >= 0,
                        "Optimist::Optimizer::Adam::decay(...): invalid input "
                        "detected.");
        this->m_decay = t_decay;
      }

      /**
       * Enable or disable the AMSGrad correction (monotone second moment).
       */
      void amsgrad_mode(const bool t_amsgrad) {
        this->m_amsgrad = t_amsgrad;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] gradient Gradient of the cost function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool solve_impl(FunctionLambda &&function,
                      GradientLambda &&gradient,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::Optimizer::Adam::solve_impl(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        const Integer n{static_cast<Integer>(x_ini.size())};
        Dense x(Utils::to_dense(x_ini)), g(n), m(Dense::Zero(n)),
            v(Dense::Zero(n)), v_hat_max(Dense::Zero(n));
        Dense x_best(x);
        Scalar f, f_best{INFTY};

        auto evaluate = [&](const Dense &z) -> bool {
          const Vector zv{Utils::from_dense<Vector>(z)};
          FirstDerivative gv;
          if (!this->evaluate_function(function, zv, f) ||
              !this->evaluate_gradient(gradient, zv, gv)) {
            return false;
          }
          g = Utils::to_dense(gv);
          return std::isfinite(f) && g.allFinite();
        };

        Scalar beta_1_t{1}, beta_2_t{1};
        std::string notes{"Initialize"};
        try {
          OPTIMIST_ASSERT(evaluate(x), CMD "evaluation failed.");
          for (this->m_iterations = 1;
               this->m_iterations < this->m_max_iterations;
               ++this->m_iterations) {
            if (f < f_best) {
              f_best = f;
              x_best = x;
            }
            const Scalar g_norm{g.norm()};
            if (this->m_verbose && (this->m_iterations % 100 == 1)) {
              this->info(g_norm, notes);
            }
            if (g_norm < this->m_tolerance) {
              this->m_converged = true;
              x_best            = x;
              break;
            }

            // Biased moment estimates and bias corrections
            m = this->m_beta_1 * m + (1 - this->m_beta_1) * g;
            v = this->m_beta_2 * v + (1 - this->m_beta_2) * g.cwiseAbs2();
            beta_1_t *= this->m_beta_1;
            beta_2_t *= this->m_beta_2;
            const Dense m_hat{m / (1 - beta_1_t)};
            Dense v_hat{v / (1 - beta_2_t)};
            if (this->m_amsgrad) {
              v_hat_max = v_hat_max.cwiseMax(v_hat);
              v_hat     = v_hat_max;
            }

            const Scalar eta{this->m_learning_rate /
                             (1 + this->m_decay * this->m_iterations)};
            x -= eta * m_hat.cwiseQuotient(
                           (v_hat.cwiseSqrt().array() + this->m_epsilon).matrix());
            OPTIMIST_ASSERT(evaluate(x), CMD "evaluation failed.");
            notes = "-";
          }
        } catch (const std::runtime_error &) {
          notes = "Budget exhausted";
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = Utils::from_dense<Vector>(x_best);
        return this->m_converged;

#undef CMD
      }

    };  // class Adam

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_ADAM_HH
