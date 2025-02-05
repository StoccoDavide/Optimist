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

#ifndef OPTIMIST_ROOTFINDER_NEWTON_HXX
#define OPTIMIST_ROOTFINDER_NEWTON_HXX

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
    * \tparam N The dimension of the nonlinear system of equations.
    */
    template <Integer N>
    class Newton : public RootFinder<N>
    {
    public:
      using Vector   = typename RootFinder<N>::Vector;
      using Matrix   = typename RootFinder<N>::Matrix;
      using Function = typename RootFinder<N>::Function;
      using Jacobian = typename RootFinder<N>::Jacobian;
      using RootFinder<N>::solve;

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
      * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0 \f$, with \f$
      * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
      * \param[in] x_ini The initialization point.
      * \param[out] x_sol The solution point.
      * \return The convergence boolean flag.
      */
      bool solve(Vector const &x_ini, Vector &x_sol) override
      {
        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        Real residuals, step_norm;
        Vector x_old, x_new, function_old, function_new, step_old, step_new;
        Matrix jacobian;

        // Set initial iteration
        x_old = x_ini;
        this->evaluate_function(x_old, function_old);

        // Algorithm iterations
        Real tolerance_residuals{this->m_tolerance};
        Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = Integer(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Store trace
          this->store_trace(x_old);

          // Calculate step
          this->evaluate_jacobian(x_old, jacobian);
          this->m_lu.compute(jacobian);
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
            bool damped{this->damp(x_old, function_old, step_old, x_new, function_new, step_new)};
            OPTIMIST_ASSERT_WARNING(damped,
              "Optimist::RootFinder::Newton::solve(...): damping failed.");
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

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_NEWTON_HXX
