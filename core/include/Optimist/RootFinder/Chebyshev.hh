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

#ifndef OPTIMIST_ROOTFINDER_CHEBYSHEV_HH
#define OPTIMIST_ROOTFINDER_CHEBYSHEV_HH

#include "Optimist/RootFinder.hh"

namespace Optimist {
  namespace RootFinder {

    /**
     * \brief Class container for the Chebyshev's method.
     *
     * \includedoc docs/markdown/RootFinder/Chebyshev.md
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    class Chebyshev : public RootFinder<Scalar, Chebyshev<Scalar>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{true};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Chebyshev solver.
       */
      Chebyshev() {}

      /**
       * Get the Chebyshev solver name.
       * \return The Chebyshev solver name.
       */
      constexpr std::string name_impl() const {
        return "Chebyshev";
      }

      /**
       * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R}
       * \rightarrow \mathbb{R} \f$.
       * \tparam FunctionLambda Function lambda type.
       * \tparam FirstDerivativeLambda First derivative lambda type.
       * \tparam SecondDerivativeLambda Second derivative lambda type.
       * \param[in] function Function lambda.
       * \param[in] first_derivative First derivative lambda.
       * \param[in] second_derivative Second derivative lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda,
                typename FirstDerivativeLambda,
                typename SecondDerivativeLambda>
      bool solve_impl(FunctionLambda &&function,
                      FirstDerivativeLambda &&first_derivative,
                      SecondDerivativeLambda &&second_derivative,
                      Scalar x_ini,
                      Scalar &x_sol) {
#define CMD "Optimist::RootFinder::Chebyshev::solve(...): "

        // Reset internal variables
        this->reset_counters();

        // Print header
        if (this->m_verbose) {
          this->header();
        }

        // Initialize variables
        bool damped, success;
        Scalar residuals, step_norm;
        Scalar x_old, x_new, function_old, function_new, step_old, step_new;
        Scalar first_derivative_old, second_derivative_old;

        // Set initial iteration
        x_old = x_ini;
        success =
            this->evaluate_function(std::forward<FunctionLambda>(function),
                                    x_old,
                                    function_old);
        OPTIMIST_ASSERT(success,
                        CMD "function evaluation failed at iteration "
                            << this->m_iterations << ".");

        // Algorithm iterations
        Scalar tolerance_residuals{this->m_tolerance};
        Scalar tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          // Evaluate derivatives
          success = this->evaluate_first_derivative(
              std::forward<FirstDerivativeLambda>(first_derivative),
              x_old,
              first_derivative_old);
          OPTIMIST_ASSERT(success,
                          CMD "first derivative evaluation failed at iteration "
                              << this->m_iterations << ".");
          success = this->evaluate_second_derivative(
              std::forward<SecondDerivativeLambda>(second_derivative),
              x_old,
              second_derivative_old);
          OPTIMIST_ASSERT(
              success,
              CMD "second derivative evaluation failed at iteration "
                  << this->m_iterations << ".");

          // Calculate step
          if (std::abs(first_derivative_old) < Chebyshev::SQRT_EPSILON) {
            OPTIMIST_WARNING(CMD
                             "close-to-singular first derivative detected.");
            first_derivative_old = (first_derivative_old > 0)
                                       ? Chebyshev::SQRT_EPSILON
                                       : -Chebyshev::SQRT_EPSILON;
          }
          if (std::abs(second_derivative_old) < Chebyshev::SQRT_EPSILON) {
            OPTIMIST_WARNING(CMD
                             "close-to-singular second derivative detected.");
            second_derivative_old = (second_derivative_old > 0)
                                        ? Chebyshev::SQRT_EPSILON
                                        : -Chebyshev::SQRT_EPSILON;
          }
          step_old = -(function_old / first_derivative_old) *
                     (1.0 + (function_old * second_derivative_old) /
                                (first_derivative_old * first_derivative_old));
          OPTIMIST_ASSERT(
              std::isfinite(step_old),
              CMD "step " << this->m_iterations << " is not finite.");

          // Check convergence
          residuals = std::abs(function_old);
          step_norm = std::abs(step_old);
          if (this->m_verbose) {
            this->info(residuals);
          }
          if (residuals < tolerance_residuals ||
              step_norm < tolerance_step_norm) {
            this->m_converged = true;
            break;
          }

          if (this->m_damped) {
            // Relax the iteration process
            damped = this->damp(std::forward<FunctionLambda>(function),
                                x_old,
                                function_old,
                                step_old,
                                x_new,
                                function_new,
                                step_new);
            OPTIMIST_ASSERT_WARNING(damped, CMD "damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            success =
                this->evaluate_function(std::forward<FunctionLambda>(function),
                                        x_new,
                                        function_new);
            OPTIMIST_ASSERT(success,
                            CMD "function evaluation failed at iteration "
                                << this->m_iterations << ".");
          }

          // Update internal variables
          x_old        = x_new;
          function_old = function_new;
          step_old     = step_new;
        }

        // Print bottom
        if (this->m_verbose) {
          this->bottom();
        }

        // Convergence data
        x_sol = x_old;
        return this->m_converged;

#undef CMD
      }

    };  // class Chebyshev

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_CHEBYSHEV_HH
