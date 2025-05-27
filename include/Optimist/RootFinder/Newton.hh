/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
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
    * \tparam Real Scalar number type.
    * \tparam N Dimension of the root-finding problem.
    */
    template <typename Real, Integer N = 1>
    class Newton : public RootFinder<Real, N, Newton<Real, N>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{false};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      using Vector = typename RootFinder<Real, N, Newton<Real, N>>::Vector;
      using Matrix = typename RootFinder<Real, N, Newton<Real, N>>::Matrix;
      using FunctionWrapper = typename RootFinder<Real, N, Newton<Real, N>>::FunctionWrapper;
      using JacobianWrapper = typename RootFinder<Real, N, Newton<Real, N>>::JacobianWrapper;
      using RootFinder<Real, N, Newton<Real, N>>::solve;

    private:
      Eigen::FullPivLU<Matrix> m_lu; /**< LU decomposition. */

    public:
      /**
      * Class constructor for the Newton solver.
      */
      Newton() {}

      /**
      * Get the Newton solver name.
      * \return The Newton solver name.
      */
      std::string name_impl() const {return "Newton";}

      /**
      * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0 \f$, with \f$
      * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
      * \param[in] function Function wrapper.
      * \param[in] jacobian Jacobian wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(FunctionWrapper function, JacobianWrapper jacobian, Vector const & x_ini, Vector & x_sol)
      {
        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped;
        Real residuals, step_norm;
        Vector x_old, x_new, function_old, function_new, step_old, step_new;
        Matrix jacobian_old;

        // Set initial iteration
        x_old = x_ini;
        this->evaluate_function(function, x_old, function_old);

        // Algorithm iterations
        Real tolerance_residuals{this->m_tolerance};
        Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = static_cast<Integer>(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Store trace
          this->store_trace(x_old);

          // Calculate step
          this->evaluate_jacobian(jacobian, x_old, jacobian_old);
          this->m_lu.compute(jacobian_old);
          OPTIMIST_ASSERT(this->m_lu.rank() == N,
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
            damped = this->damp(function, x_old, function_old, step_old, x_new, function_new, step_new);
            OPTIMIST_ASSERT_WARNING(damped,
              "Optimist::RootFinder::Newton::solve(...): damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            this->evaluate_function(function, x_new, function_new);
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
      }

    }; // class Newton

    /**
    * \brief Class container for the Newton's method.
    *
    * \includedoc docs/markdown/RootFinder/NewtonScalar.md
    *
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Newton<Real, 1> : public RootFinder<Real, 1, Newton<Real, 1>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{false};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      // Function types
      using FunctionWrapper         = typename RootFinder<Real, 1, Newton<Real, 1>>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename RootFinder<Real, 1, Newton<Real, 1>>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename RootFinder<Real, 1, Newton<Real, 1>>::SecondDerivativeWrapper;

      /**
      * Class constructor for the Newton solver.
      */
      Newton<Real, 1>() {}

      /**
      * Get the Newton solver name.
      * \return The Newton solver name.
      */
      std::string name_impl() const {return "Newton";}

      /**
      * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R} \rightarrow \mathbb{R} \f$.
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(FunctionWrapper function, FirstDerivativeWrapper first_derivative, Real x_ini,
        Real & x_sol)
      {
        #define CMD "Optimist::RootFinder::Newton::solve(...): "

        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped;
        Real residuals, step_norm;
        Real x_old, x_new, function_old, function_new, step_old, step_new;
        Real first_derivative_old;

        // Set initial iteration
        x_old = x_ini;
        this->evaluate_function(function, x_old, function_old);

        // Algorithm iterations
        Real tolerance_residuals{this->m_tolerance};
        Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = static_cast<Integer>(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Store trace
          this->store_trace(x_old);

          // Evaluate first derivative
          this->evaluate_first_derivative(first_derivative, x_old, first_derivative_old);

          // Calculate step
          if (std::abs(first_derivative_old) < EPSILON_LOW) {
            OPTIMIST_WARNING( CMD "singular first derivative detected.");
            first_derivative_old = (first_derivative_old > 0.0) ? EPSILON_LOW : -EPSILON_LOW;
          }
          step_old = -function_old/first_derivative_old;
          OPTIMIST_ASSERT(std::isfinite(step_old), CMD "step " << this->m_iterations << " is not finite.");

          // Check convergence
          residuals = std::abs(function_old);
          step_norm = std::abs(step_old);
          if (this->m_verbose) {this->info(residuals);}
          if (residuals < tolerance_residuals || step_norm < tolerance_step_norm) {
            this->m_converged = true;
            break;
          }

          if (this->m_damped) {
            // Relax the iteration process
            damped = this->damp(function, x_old, function_old, step_old, x_new, function_new, step_new);
            OPTIMIST_ASSERT_WARNING(damped, CMD "damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            this->evaluate_function(function, x_new, function_new);
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
