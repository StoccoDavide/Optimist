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

#ifndef OPTIMIST_SCALAR_ROOT_FINDER_CHEBYSHEV_HXX
#define OPTIMIST_SCALAR_ROOT_FINDER_CHEBYSHEV_HXX

namespace Optimist
{
  namespace ScalarRootFinder
  {

    /*\
     |    ____ _          _               _
     |   / ___| |__   ___| |__  _   _ ___| |__   _____   __
     |  | |   | '_ \ / _ \ '_ \| | | / __| '_ \ / _ \ \ / /
     |  | |___| | | |  __/ |_) | |_| \__ \ | | |  __/\ V /
     |   \____|_| |_|\___|_.__/ \__, |___/_| |_|\___| \_/
     |                          |___/
    \*/

    /**
    * \brief Class container for the Chebyshev's method.
    *
    * \includedoc docs/markdown/ScalarRootFinder/Chebyshev.md
    *
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Chebyshev : public ScalarRootFinder<Real, Chebyshev<Real>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{true};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      // Function types
      using FunctionWrapper         = typename ScalarRootFinder<Real, Chebyshev<Real>>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename ScalarRootFinder<Real, Chebyshev<Real>>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename ScalarRootFinder<Real, Chebyshev<Real>>::SecondDerivativeWrapper;

      /**
      * Class constructor for the Chebyshev solver.
      */
      Chebyshev() {}

      /**
      * Get the Chebyshev solver name.
      * \return The Chebyshev solver name.
      */
      std::string name_impl() const {return "Chebyshev";}

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
        #define CMD "Optimist::ScalarRootFinder::Chebyshev::solve(...): "

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
        for (this->m_iterations = Integer(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
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
          step_old = -(function_old / first_derivative_old) * (1.0 + (function_old * second_derivative_old) /
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

    }; // class Chebyshev

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_CHEBYSHEV_HXX
