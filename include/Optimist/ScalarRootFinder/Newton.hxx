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

#ifndef OPTIMIST_SCALAR_ROOT_FINDER_NEWTON_HXX
#define OPTIMIST_SCALAR_ROOT_FINDER_NEWTON_HXX

namespace Optimist
{
  namespace ScalarRootFinder
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
    * \includedoc docs/markdown/ScalarRootFinder/Newton.md
    */
    class Newton : public ScalarRootFinder
    {
    public:
      // Function types
      using Function         = typename ScalarRootFinder::Function;        /**< Function type. */
      using FirstDerivative  = typename ScalarRootFinder::FirstDerivative; /**<  First derivative type. */
      using SecondDerivative = typename ScalarRootFinder::SecondDerivative; /**< Second derivative type. */

      /**
      * Class constructor for the Newton solver.
      */
      Newton() {}

      /**
      * Get the Newton solver name.
      * \return The Newton solver name.
      */
      std::string name() const override {return "Newton";}

      /**
      * Check if the Broyden solver is able to solve the problem with the given input.
      * \return The check boolean flag.
      */
      bool check() const override
      {
        if (this->m_function != nullptr && this->m_first_derivative != nullptr) {
          return true;
        } else {
          return false;
        }
      }

      /**
      * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R} \rightarrow \mathbb{R} \f$.
      * \param[in] x_ini The initialization point.
      * \param[out] x_sol The solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Real const &x_ini, Real &x_sol) override
      {
        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        Real residuals, step_norm;
        Real x_old, x_new, function_old, function_new, step_old, step_new;
        Real first_derivative;

        // Set initial iteration
        x_old = x_ini;
        this->evaluate_function(x_old, function_old);

        // Algorithm iterations
        Real tolerance_residuals{this->m_tolerance};
        Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = Integer(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Store trace
          this->evaluate_first_derivative(x_new, first_derivative);
          this->store_trace(x_old);

          // Calculate step
          OPTIMIST_ASSERT(std::abs(first_derivative) > EPSILON_LOW,
            "Optimist::ScalarRootFinder::Newton::solve(...): singular first derivative detected.");
          step_old = -function_old/first_derivative;

          // Check convergence
          residuals = std::abs(function_old);
          step_norm = std::abs(step_old);
          if (this->m_verbose) {this->info(residuals);}
          if (residuals < tolerance_residuals || step_norm < tolerance_step_norm) {
            this->m_converged = true;
            break;
          }

          if (this->m_damped)
          {
            // Relax the iteration process
            bool damped{this->damp(x_old, function_old, step_old, x_new, function_new, step_new)};
            OPTIMIST_ASSERT_WARNING(damped,
              "Optimist::ScalarRootFinder::Newton::solve(...): damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            this->evaluate_function(x_new, function_new);
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

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_NEWTON_HXX
