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

#ifndef OPTIMIST_ROOTFINDER_QUASINEWTON_HH
#define OPTIMIST_ROOTFINDER_QUASINEWTON_HH

#include "Optimist/RootFinder.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |    ___                  _ _   _               _
     |   / _ \ _   _  __ _ ___(_) \ | | _____      _| |_ ___  _ __
     |  | | | | | | |/ _` / __| |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
     |  | |_| | |_| | (_| \__ \ | |\  |  __/\ V  V /| || (_) | | | |
     |   \__\_\\__,_|\__,_|___/_|_| \_|\___| \_/\_/  \__\___/|_| |_|
     |
    \*/

    /**
     * \brief Class container for the QuasiNewton's method.
     *
     * \includedoc docs/markdown/RootFinder/QuasiNewton.md
     *
     * \tparam Real Scalar number type.
     * \tparam N Dimension of the root-finding problem.
     */
    template <typename Real, Integer N, typename DerivedSolver>
    class QuasiNewton : public RootFinder<Real, N, DerivedSolver, true>
    {
    public:
      OPTIMIST_BASIC_CONSTANTS(Real)

      using typename RootFinder<Real, N, DerivedSolver, true>::Vector;
      using typename RootFinder<Real, N, DerivedSolver, true>::Matrix;

      /**
       * QuasiNewton solver type enumeration.
       */
      using Method = enum class Method : Integer {GOOD = 0, BAD = 1, COMBINED = 2};

    private:
      Method m_method{Method::COMBINED}; /**< QuasiNewton solver type. */

    public:
      /**
       * Class constructor for the QuasiNewton solver.
       */
      QuasiNewton() {}

      /**
       * Get the QuasiNewton solver name.
       * \return The QuasiNewton solver name.
       */
      std::string name_impl() const
      {
        return static_cast<const DerivedSolver *>(this)->name_impl();
      }

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
        #define CMD "Optimist::RootFinder::QuasiNewton::solve(...): "

        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped, success;
        Real residuals, step_norm;
        Vector x_old, x_new, function_old, function_new, step_old, step_new, delta_x_old, delta_x_new,
          delta_function_old, delta_function_new;
        Matrix jacobian_old, jacobian_new;

        // Set initial iteration
        x_old = x_ini;
        success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old, function_old);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the initial point.");
        success = this->evaluate_jacobian(jacobian, x_old, jacobian_old);
        OPTIMIST_ASSERT(success,
          CMD "jacobian evaluation failed at the initial point.");

        // Algorithm iterations
        Real tolerance_residuals{this->m_tolerance};
        Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1; this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Calculate step
          step_old = -jacobian_old * function_old;

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
              "Optimist::RootFinder::QuasiNewton::solve(...): damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            success = this->evaluate_function(std::forward<FunctionLambda>(function), x_new, function_new);
            OPTIMIST_ASSERT(success,
              CMD "function evaluation failed at iteration " << this->m_iterations << ".");
          }

          // Update jacobian approximation
          delta_x_new = x_new - x_old;
          delta_function_new = function_new - function_old;
          this->update(
            delta_x_old, delta_function_old, jacobian_old, // Old step data
            delta_x_new, delta_function_new, function_new, jacobian_new  // New step data
          );

          // Update internal variables
          x_old              = x_new;
          function_old       = function_new;
          delta_x_old        = delta_x_new;
          delta_function_old = delta_function_new;
          step_old           = step_new;
          jacobian_old       = jacobian_new;
        }

        // Print bottom
        if (this->m_verbose) {this->bottom();}

        // Convergence data
        x_sol = x_old;
        return this->m_converged;

        #undef CMD
      }

      /**
       * Jacobian approximation update rule for the QuasiNewton's method.
       * \param[in] delta_x_old Old difference between points.
       * \param[in] delta_function_old Old difference between function values.
       * \param[in] jacobian_old Old jacobian approximation.
       * \param[in] delta_x_new New difference between points.
       * \param[in] delta_function_new New difference between function values.
       * \param[in] function_new New function value.
       * \param[out] jacobian_new New jacobian approximation.
       */
      void update(
        Vector const & delta_x_old, Vector const & delta_function_old, Matrix const & jacobian_old,
        Vector const & delta_x_new, Vector const & delta_function_new, Vector const & function_new,
        Matrix       & jacobian_new
      ) {
        static_cast<DerivedSolver *>(this)->update_impl(
          delta_x_old, delta_function_old, jacobian_old, // Old step data
          delta_x_new, delta_function_new, function_new, jacobian_new  // New step data
        );
      }

    }; // class QuasiNewton

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_QUASINEWTON_HH
