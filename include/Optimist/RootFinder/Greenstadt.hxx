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

#ifndef OPTIMIST_ROOTFINDER_GREENSTADT_HXX
#define OPTIMIST_ROOTFINDER_GREENSTADT_HXX

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |    ____                         _            _ _
     |   / ___|_ __ ___  ___ _ __  ___| |_ __ _  __| | |_
     |  | |  _| '__/ _ \/ _ \ '_ \/ __| __/ _` |/ _` | __|
     |  | |_| | | |  __/  __/ | | \__ \ || (_| | (_| | |_
     |   \____|_|  \___|\___|_| |_|___/\__\__,_|\__,_|\__|
     |
    \*/

    /**
    * \brief Class container for the Greenstadt's method.
    *
    * \includedoc docs/markdown/RootFinder/Greenstadt.md
    *
    * \tparam N The dimension of the nonlinear system of equations.
    */
    template <Integer N>
    class Greenstadt : public RootFinder<N>
    {
    public:
      using Method = enum class Method : Integer {ONE = 1, TWO = 2}; /**< Greenstadt solver type. */
      using Vector   = typename RootFinder<N>::Vector;
      using Matrix   = typename RootFinder<N>::Matrix;
      using Function = typename RootFinder<N>::Function;
      using Jacobian = typename RootFinder<N>::Jacobian;
      using RootFinder<N>::solve;

    private:
      Method m_method{Method::ONE}; /**< Greenstadt solver type. */

    public:
      /**
      * Class constructor for the Greenstadt solver.
      */
      Greenstadt() {}

      /**
      * Get the Greenstadt solver name.
      * \return The Greenstadt solver name.
      */
      std::string name() const override
      {
        std::ostringstream os;
        os << "Greenstadt";
        if (this->m_method == Method::ONE) {
          os << "1";
        } else if (this->m_method == Method::TWO) {
          os << "2";
        }
        return os.str();
      }

      /**
      * Get the enumeration type of the Greenstadt solver method.
      * \return The Greenstadt solver enumeration type.
      */
      Method method() const {return this->m_method;}

      /**
      * Set the enumeration type of the Greenstadt solver method.
      * \param[in] t_method The Greenstadt solver method enumeration type.
      */
      void method(Method t_method) {this->m_method = t_method;}

      /**
      * Enable the \em Greenstadt1 solver method.
      */
      void enable_one_method() {this->m_method = Method::ONE;}

      /**
      * Enable the \em Greenstadt2 solver method.
      */
      void enable_two_method() {this->m_method = Method::TWO;}

      /**
      * Set the Greenstadt solver method.
      * \param[in] t_method The Greenstadt solver method enumeration type.
      */
      void set_method(Method t_method) {this->m_method = t_method;}

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
        Vector x_old, x_new, function_old, function_new, step_old, step_new, delta_x_old, delta_x_new,
          delta_function_old, delta_function_new;
        Matrix jacobian_old, jacobian_new;

        // Set initial iteration
        x_old = x_ini;
        this->evaluate_function(x_old, function_old);
        this->evaluate_jacobian(x_old, jacobian_old);

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
            bool damped{this->damp(x_old, function_old, step_old, x_new, function_new, step_new)};
            OPTIMIST_ASSERT_WARNING(damped,
              "Optimist::RootFinder::Greenstadt::solve(...): damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            this->evaluate_function(x_new, function_new);
          }

          // Update jacobian approximation
          delta_x_new = x_new - x_old;
          delta_function_new = function_new - function_old;
          this->update(
            jacobian_old, // Old step data
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
      }

      /**
      * Jacobian approximation update rule for the Greenstadt's method.
      * \param[in] jacobian_old Old jacobian approximation.
      * \param[in] delta_x_new New difference between points.
      * \param[in] delta_function_new New difference between function values.
      * \param[out] function_new New function value.
      * \param[out] jacobian_new New jacobian approximation.
      */
      void update(
        Matrix const &jacobian_old,  Vector const &delta_x_new, Vector const &delta_function_new,
        Vector const &function_new,  Matrix &jacobian_new
      ) {
        if (this->m_method == Method::ONE) {
          // Greenstadt's 1st method
          // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C = F1;
          jacobian_new = jacobian_old - (jacobian_old*delta_function_new-delta_x_new)/(function_new.transpose()*delta_function_new)*function_new.transpose();
        } else if (this->m_method == Method::TWO) {
          // Greenstadt's 2nd method
          // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C  = J0'*J0*DF1;
          Vector C(jacobian_old.transpose()*jacobian_old*delta_function_new);
          jacobian_new = jacobian_old - (jacobian_old*delta_function_new-delta_x_new)/(C.transpose()*delta_function_new)*C.transpose();
        }
      }

    }; // class Greenstadt

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_GREENSTADT_HXX
