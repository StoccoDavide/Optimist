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

#ifndef OPTIMIST_ROOTFINDER_BROYDEN_HXX
#define OPTIMIST_ROOTFINDER_BROYDEN_HXX

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |   ____                      _
     |  | __ ) _ __ ___  _   _  __| | ___ _ __
     |  |  _ \| '__/ _ \| | | |/ _` |/ _ \ '_ \
     |  | |_) | | | (_) | |_| | (_| |  __/ | | |
     |  |____/|_|  \___/ \__, |\__,_|\___|_| |_|
     |                   |___/
    \*/

    /**
    * \brief Class container for the Broyden's method.
    *
    * \includedoc docs/markdown/RootFinder/Broyden.md
    *
    * \tparam N Dimension of the root-finding problem.
    */
    template <Integer N>
    class Broyden : public RootFinder<N, Broyden<N>>
    {
    static constexpr bool requires_function          = true;
    static constexpr bool requires_first_derivative  = false;
    static constexpr bool requires_second_derivative = false;

    public:
      using Method = enum class Method : Integer {GOOD = 0, BAD = 1, COMBINED = 2}; /**< Broyden solver type. */
      using Vector   = typename RootFinder<N, Broyden<N>>::Vector;
      using Matrix   = typename RootFinder<N, Broyden<N>>::Matrix;
      using Function = typename RootFinder<N, Broyden<N>>::Function;
      using Jacobian = typename RootFinder<N, Broyden<N>>::Jacobian;
      using RootFinder<N, Broyden<N>>::solve;

    private:
      Method m_method{Method::COMBINED}; /**< Broyden solver type. */

    public:
      /**
      * Class constructor for the Broyden solver.
      */
      Broyden() {}

      /**
      * Get the Broyden solver name.
      * \return The Broyden solver name.
      */
      std::string name_impl() const
      {
        std::ostringstream os;
        os << "Broyden";
        if (this->m_method == Method::GOOD) {
          os << "Good";
        } else if (this->m_method == Method::BAD) {
          os << "Bad";
        } else if (this->m_method == Method::COMBINED) {
          os << "Combined";
        }
        return os.str();
      }

      /**
      * Get the enumeration type of the Broyden solver method.
      * \return The Broyden solver enumeration type.
      */
      Method method() const {return this->m_method;}

      /**
      * Set the enumeration type of the Broyden solver method.
      * \param[in] t_method The Broyden solver enumeration type.
      */
      void method(Method t_method) {this->m_method = t_method;}

      /**
      * Enable the \em good Broyden solver method.
      */
      void enable_good_method() {this->m_method = Method::GOOD;}

      /**
      * Enable the \em bad Broyden solver method.
      */
      void enable_bad_method() {this->m_method = Method::BAD;}

      /**
      * Enable the \em combined Broyden solver method.
      */
      void enable_combined_method() {this->m_method = Method::COMBINED;}

      /**
      * Set the Broyden solver type.
      * \param[in] t_method The Broyden solver type enumeration.
      */
      void set_method(Method t_method) {this->m_method = t_method;}

      /**
      * Solve the nonlinear system of equations \f$ \mathbf{f}(\mathbf{x}) = 0 \f$, with \f$
      * \mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^n \f$.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(Function t_function, Jacobian t_jacobian, Vector const &x_ini, Vector &x_sol)
      {
        // Setup internal variables
        this->reset();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped;
        Real residuals, step_norm;
        Vector x_old, x_new, function_old, function_new, step_old, step_new, delta_x_old, delta_x_new,
          delta_function_old, delta_function_new;
        Matrix jacobian_old, jacobian_new;

        // Set initial iteration
        x_old = x_ini;
        this->evaluate_function(t_function, x_old, function_old);
        this->evaluate_jacobian(t_jacobian, x_old, jacobian_old);

        // Algorithm iterations
        Real tolerance_residuals{this->m_tolerance};
        Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = Integer(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Store trace
          this->store_trace(x_old);

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
            damped = this->damp(t_function, x_old, function_old, step_old, x_new, function_new, step_new);
            OPTIMIST_ASSERT_WARNING(damped,
              "Optimist::RootFinder::Broyden::solve(...): damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            this->evaluate_function(t_function, x_new, function_new);
          }

          // Update jacobian approximation
          delta_x_new = x_new - x_old;
          delta_function_new = function_new - function_old;
          this->update(
            delta_x_old, delta_function_old, jacobian_old, // Old step data
            delta_x_new, delta_function_new, jacobian_new  // New step data
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
      }

      /**
      * Jacobian approximation update rule for the Broyden's method.
      * \param[in] delta_x_old Old difference between points.
      * \param[in] delta_function_old Old difference between function values.
      * \param[in] jacobian_old Old jacobian approximation.
      * \param[in] delta_x_new New difference between points.
      * \param[in] delta_function_new New difference between function values.
      * \param[out] jacobian_new New jacobian approximation.
      */
      void update(
        Vector const &delta_x_old, Vector const &delta_function_old, Matrix const &jacobian_old,
        Vector const &delta_x_new, Vector const &delta_function_new, Matrix       &jacobian_new
      ) {
        Vector tmp_1(jacobian_old * delta_function_new);
        Real tmp_2{delta_function_new.squaredNorm()};
        // Selection criteria: |(dx_new'*dx_old) / (dx_new'*J_old*dF_new)| < |(dF_new'*dF_old) / (dF_new'*dF_new)|
        if (this->m_method == Method::COMBINED || this->m_method == Method::GOOD || this->iterations() < Integer(2) ||
            std::abs(delta_x_new.transpose() * delta_x_old) / std::abs(delta_x_new.transpose() * tmp_1)
            < std::abs(delta_function_new.transpose() * delta_function_old) / tmp_2) {
          // Broyden's Good solver: J_new = J_old - (J_old*dF_new-dx_new) / (C'*dF_new)*C', where C = J_old'*dx_new;
          Vector C_g(jacobian_old.transpose() * delta_x_new);
          jacobian_new = jacobian_old - (tmp_1 - delta_x_new) / (C_g.transpose() * delta_function_new) * C_g.transpose();
        } else {
          // Broyden's Bad solver: J_new = J_old - (J_old*dF_old-dx_new) / (C'*dF_old)*C', where C = dF_old;
          jacobian_new = jacobian_old - (tmp_1 - delta_x_new) / tmp_2 * delta_function_old.transpose();
        }
      }

    }; // class Broyden

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_BROYDEN_HXX
