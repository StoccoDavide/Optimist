/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef OPTIMIST_ROOTFINDER_NEWTONRAPHSON_HH
#define OPTIMIST_ROOTFINDER_NEWTONRAPHSON_HH

#include "Optimist/RootFinder.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |   _   _               _              ____             _
     |  | \ | | _____      _| |_ ___  _ __ |  _ \ __ _ _ __ | |__  ___  ___  _ __
     |  |  \| |/ _ \ \ /\ / / __/ _ \| '_ \| |_) / _` | '_ \| '_ \/ __|/ _ \| '_ \
     |  | |\  |  __/\ V  V /| || (_) | | | |  _ < (_| | |_) | | | \__ \ (_) | | | |
     |  |_| \_|\___| \_/\_/  \__\___/|_| |_|_| \_\__,_| .__/|_| |_|___/\___/|_| |_|
     |                                                |_|
    \*/

    /**
    * \brief Class container for the scalar Newton-Raphson method.
    *
    * \includedoc docs/markdown/RootFinder/NewtonRaphson.md
    *
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class NewtonRaphson : public RootFinder<Real, 1, NewtonRaphson<Real>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{false};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      // Function types
      using FunctionWrapper         = typename RootFinder<Real, 1, NewtonRaphson<Real>>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename RootFinder<Real, 1, NewtonRaphson<Real>>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename RootFinder<Real, 1, NewtonRaphson<Real>>::SecondDerivativeWrapper;

      /**
      * Class constructor for the Newton solver.
      */
      NewtonRaphson() {}

      /**
      * Get the Newton solver name.
      * \return The Newton solver name.
      */
      std::string name_impl() const {return "NewtonRaphson";}

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
        #define CMD "Optimist::RootFinder::NewtonRaphson::solve(...): "

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

    }; // class NewtonRaphson

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_NEWTONRAPHSON_HH
