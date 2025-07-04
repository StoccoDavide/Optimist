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

#ifndef OPTIMIST_ROOTFINDER_HALLEY_HH
#define OPTIMIST_ROOTFINDER_HALLEY_HH

#include "Optimist/RootFinder.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |   _   _       _ _
     |  | | | | __ _| | | ___ _   _
     |  | |_| |/ _` | | |/ _ \ | | |
     |  |  _  | (_| | | |  __/ |_| |
     |  |_| |_|\__,_|_|_|\___|\__, |
     |                        |___/
    \*/

    /**
    * \brief Class container for the Halley's method.
    *
    * \includedoc docs/markdown/RootFinder/Halley.md
    *
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Halley : public RootFinder<Real, 1, Halley<Real>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{true};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      // Function types
      using FunctionWrapper         = typename RootFinder<Real, 1, Halley<Real>>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename RootFinder<Real, 1, Halley<Real>>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename RootFinder<Real, 1, Halley<Real>>::SecondDerivativeWrapper;

      /**
      * Class constructor for the Halley solver.
      */
      Halley() {}

      /**
      * Get the Halley solver name.
      * \return The Halley solver name.
      */
      std::string name_impl() const {return "Halley";}

      /**
      * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R} \rightarrow \mathbb{R} \f$.
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] second_derivative Second derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(FunctionWrapper function, FirstDerivativeWrapper first_derivative,
        SecondDerivativeWrapper second_derivative, Real x_ini, Real & x_sol)
      {
        #define CMD "Optimist::RootFinder::Halley::solve(...): "

        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped;
        Real residuals, step_norm;
        Real x_old, x_new, function_old, function_new, step_old, step_new;
        Real first_derivative_old, second_derivative_old;

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

          // Evaluate derivatives
          this->evaluate_first_derivative(first_derivative, x_old, first_derivative_old);
          this->evaluate_second_derivative(second_derivative, x_old, second_derivative_old);

          // Calculate step
          if (std::abs(first_derivative_old) < EPSILON_LOW) {
            OPTIMIST_WARNING( CMD "singular first derivative detected.");
            first_derivative_old = (first_derivative_old > 0.0) ? EPSILON_LOW : -EPSILON_LOW;
          }
          if (std::abs(second_derivative_old) < EPSILON_LOW) {
            OPTIMIST_WARNING( CMD "singular second derivative detected.");
            second_derivative_old = (second_derivative_old > 0.0) ? EPSILON_LOW : -EPSILON_LOW;
          }
          step_old = -(function_old / first_derivative_old) * (1.0 - (function_old * second_derivative_old) /
            (first_derivative_old * first_derivative_old));
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

    }; // class Halley

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_HALLEY_HH
