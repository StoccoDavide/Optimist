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

#ifndef OPTIMIST_ROOTFINDER_NEWTON_HH
#define OPTIMIST_ROOTFINDER_NEWTON_HH

#include "Optimist/RootFinder.hh"

namespace Optimist {
  namespace RootFinder {

    /**
     * \brief Class container for the Newton's method.
     *
     * \includedoc docs/markdown/RootFinder/Newton.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class Newton : public RootFinder<Vector, Newton<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename TypeTrait<Vector>::Scalar;
      using typename RootFinder<Vector, Newton<Vector>>::FirstDerivative;
      using Factorization =
          std::conditional_t<VectorTrait::IsSparse,
                             Eigen::SparseLU<FirstDerivative>,
                             Eigen::FullPivLU<FirstDerivative>>;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Factorization m_lu; /**< LU decomposition. */

     public:
      /**
       * Class constructor for the Newton solver.
       */
      Newton() {}

      /**
       * Get the Newton solver name.
       * \return The Newton solver name.
       */
      constexpr std::string name_impl() const {
        return "Newton";
      }

      /**
       * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0
       * \f$, with \f$
       * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
       * \tparam FunctionLambda Function lambda type.
       * \tparam JacobianLambda Jacobian lambda type.
       * \param[in] function Function lambda.
       * \param[in] jacobian Jacobian lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename JacobianLambda>
      bool solve_impl(FunctionLambda &&function,
                      JacobianLambda &&jacobian,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::RootFinder::Newton::solve(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        Vector x(x_ini), f, step, x_new, f_new, step_new;
        FirstDerivative J;
        OPTIMIST_ASSERT(this->evaluate_function(function, x, f),
                        CMD "function evaluation failed at the initial point.");

        const Scalar tolerance_residuals{this->m_tolerance};
        const Scalar tolerance_step{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          // Newton step: J h = -f
          OPTIMIST_ASSERT(this->evaluate_jacobian(jacobian, x, J),
                          CMD "jacobian evaluation failed at iteration "
                              << this->m_iterations << ".");
          this->m_lu.compute(J);
          if constexpr (VectorTrait::IsSparse) {
            OPTIMIST_ASSERT(this->m_lu.info() == Eigen::Success,
                            CMD "jacobian factorization failed.");
          } else {
            OPTIMIST_ASSERT(this->m_lu.isInvertible(),
                            CMD "singular jacobian detected.");
          }
          step = this->m_lu.solve(f);
          step = -step;

          // Convergence
          if (this->m_verbose) {
            this->info(f.norm());
          }
          if (f.norm() < tolerance_residuals || step.norm() < tolerance_step) {
            this->m_converged = true;
            break;
          }

          // Update
          if (this->m_damped) {
            this->damp(function, x, f, step, x_new, f_new, step_new);
          } else {
            x_new = x + step;
            OPTIMIST_ASSERT(this->evaluate_function(function, x_new, f_new),
                            CMD "function evaluation failed at iteration "
                                << this->m_iterations << ".");
          }
          x = x_new;
          f = f_new;
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = x;
        return this->m_converged;

#undef CMD
      }

    };  // class Newton

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_NEWTON_HH
