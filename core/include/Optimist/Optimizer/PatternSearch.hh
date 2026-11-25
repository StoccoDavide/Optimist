/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_OPTIMIZER_PATTERN_SEARCH_HH
#define OPTIMIST_OPTIMIZER_PATTERN_SEARCH_HH

#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the Hooke-Jeeves pattern search.
     *
     * \includedoc docs/markdown/Optimizer/PatternSearch.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class PatternSearch : public Optimizer<Vector, PatternSearch<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename Vector::Scalar;
      using Dense       = Utils::DenseVector<Vector>;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_step{1.0}; /**< Initial stencil step. */
      Scalar m_rho{0.5};  /**< Stencil step reduction factor in (0, 1). */

     public:
      /**
       * Class constructor for the pattern search solver.
       */
      PatternSearch() {}

      /**
       * Get the pattern search solver name.
       */
      constexpr std::string name_impl() const {
        return "PatternSearch";
      }

      /**
       * Get the initial stencil step.
       */
      Scalar step() const {
        return this->m_step;
      }

      /**
       * Set the initial stencil step.
       */
      void step(const Scalar t_step) {
        OPTIMIST_ASSERT(t_step > 0,
                        "Optimist::Optimizer::PatternSearch::step(...): "
                        "invalid input detected.");
        this->m_step = t_step;
      }

      /**
       * Get the stencil step reduction factor.
       */
      Scalar rho() const {
        return this->m_rho;
      }

      /**
       * Set the stencil step reduction factor \f$ 0 < \rho < 1 \f$.
       */
      void rho(const Scalar t_rho) {
        OPTIMIST_ASSERT(t_rho > 0 && t_rho < 1,
                        "Optimist::Optimizer::PatternSearch::rho(...): "
                        "invalid input detected.");
        this->m_rho = t_rho;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] x_ini Initialization point.
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
        Dense base(Utils::to_dense(x_ini));
        Dense sign(Dense::Ones(n));  // Preferred direction per coordinate
        Scalar h{this->m_step};
        Scalar f_base{this->eval(function, base)};
        std::string notes{"Initialize"};

        try {
          for (this->m_iterations = 1;
               this->m_iterations < this->m_max_iterations;
               ++this->m_iterations) {
            if (this->m_verbose) {
              this->info(h, notes);
            }
            if (h < this->m_tolerance * (1 + base.template lpNorm<Eigen::Infinity>())) {
              this->m_converged = true;
              break;
            }

            // Exploratory moves around the base point
            Dense x(base);
            Scalar f{f_base};
            this->explore(function, x, f, sign, h);
            if (!(f < f_base)) {
              h *= this->m_rho;
              notes = "Shrink";
              continue;
            }

            // Pattern moves: keep jumping along the successful direction as
            // long as the exploration around the jump improves
            notes = "Explore";
            while (true) {
              Dense jump{2 * x - base};
              base   = x;
              f_base = f;
              Scalar f_jump{this->eval(function, jump)};
              this->explore(function, jump, f_jump, sign, h);
              if (f_jump < f_base) {
                x     = jump;
                f     = f_jump;
                notes = "Pattern";
              } else {
                break;
              }
            }
          }
        } catch (const std::runtime_error &) {
          notes = "Budget exhausted";
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = Utils::from_dense<Vector>(base);
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
       * Coordinate-wise exploration with step h, trying the preferred
       * direction first.
       */
      template <typename FunctionLambda>
      void explore(FunctionLambda &function,
                   Dense &x,
                   Scalar &f,
                   Dense &sign,
                   const Scalar h) {
        for (Integer j{0}; j < static_cast<Integer>(x.size()); ++j) {
          const Scalar x_j{x(j)};
          x(j) = x_j + sign(j) * h;
          Scalar f_try{this->eval(function, x)};
          if (f_try < f) {
            f = f_try;
            continue;
          }
          x(j)  = x_j - sign(j) * h;
          f_try = this->eval(function, x);
          if (f_try < f) {
            f       = f_try;
            sign(j) = -sign(j);
            continue;
          }
          x(j) = x_j;
        }
      }

    };  // class PatternSearch

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_PATTERN_SEARCH_HH
