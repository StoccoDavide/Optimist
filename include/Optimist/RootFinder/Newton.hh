/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_ROOTFINDER_NEWTON_HH
#define OPTIMIST_ROOTFINDER_NEWTON_HH

#include "Optimist/RootFinder.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |   _   _               _
     |  | \ | | _____      _| |_ ___  _ __
     |  |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
     |  | |\  |  __/\ V  V /| || (_) | | | |
     |  |_| \_|\___| \_/\_/  \__\___/|_| |_|
     |
    \*/

    /**
     * \brief Class container for the Newton's method.
     *
     * \includedoc docs/markdown/RootFinder/Newton.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
    requires TypeTrait<Vector>::IsEigen
    class Newton : public RootFinder<Vector, Vector, Newton<Vector>>
    {
    public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar = typename Vector::Scalar;
      using typename RootFinder<Vector, Vector, Newton<Vector>>::Vector;
      using typename RootFinder<Vector, Vector, Newton<Vector>>::Matrix;
      using Factorization = std::conditional_t<VectorTrait::IsSparse, Eigen::SparseLU<Matrix>,
        Eigen::FullPivLU<Matrix>>;

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
      constexpr std::string name_impl() const {return "Newton";}

      /**
       * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0 \f$, with \f$
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
      bool solve_impl(FunctionLambda && function, JacobianLambda && jacobian, Vector const & x_ini,
        Vector & x_sol)
      {
        #define CMD "Optimist::RootFinder::Newton::solve(...): "

        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped, success;
        Scalar residuals, step_norm;
        Vector x_old, x_new, function_old, function_new, step_old, step_new;
        Matrix jacobian_old;

        // Set initial iteration
        x_old = x_ini;
        success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old, function_old);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the initial point.");

        // Algorithm iterations
        Scalar tolerance_residuals{this->m_tolerance};
        Scalar tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1; this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Calculate step
          success = this->evaluate_jacobian(jacobian, x_old, jacobian_old);
          OPTIMIST_ASSERT(success,
            CMD "jacobian evaluation failed at iteration " << this->m_iterations << ".");
          this->m_lu.compute(jacobian_old);
          OPTIMIST_ASSERT(this->m_lu.rank() == x_ini.size(),
            "Optimist::RootFinder::Newton::solve(...): singular Jacobian detected.");
          step_old = -this->m_lu.solve(function_old);

          // Check convergence
          residuals = function_old.norm();
          step_norm = step_old.norm();
          if (this->m_verbose) {this->info(residuals);}
          if (residuals < tolerance_residuals || step_norm < tolerance_step_norm) {
            this->m_converged = true;
            break;
          }

          if (this->m_damped)
          {
            // Relax the iteration process
            damped = this->damp(std::forward<FunctionLambda>(function), x_old, function_old, step_old, x_new, function_new, step_new);
            OPTIMIST_ASSERT_WARNING(damped,
              "Optimist::RootFinder::Newton::solve(...): damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            success = this->evaluate_function(std::forward<FunctionLambda>(function), x_new, function_new);
            OPTIMIST_ASSERT(success,
              CMD "function evaluation failed at iteration " << this->m_iterations << ".");
          }

          // Update internal variables
          x_old        = x_new;
          function_old = function_new;
          step_old     = step_new;
        }

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        x_sol = x_old;
        return this->m_converged;

        #undef CMD
      }

    }; // class Newton

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_NEWTON_HH
