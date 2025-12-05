/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco.                                                            *
 *                                                                                               *
 * The Optimist project is distributed under the BSD 2-Clause License.                           *
 *                                                                                               *
 * Davide Stocco                                                                                 *
 * University of Trento                                                                          *
 * davide.stocco@unitn.it                                                                        *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_ROOTFINDER_VARONA4_HH
#define OPTIMIST_ROOTFINDER_VARONA4_HH

#include "Optimist/RootFinder.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |  __     __                          _  _
     |  \ \   / /_ _ _ __ ___  _ __   __ _| || |
     |   \ \ / / _` | '__/ _ \| '_ \ / _` | || |_
     |    \ V / (_| | | | (_) | | | | (_| |__   _|
     |     \_/ \__,_|_|  \___/|_| |_|\__,_|  |_|
     |
    \*/

    /**
     * \brief Class container for the Varona's methods.
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    class Varona : public RootFinder<Scalar, Scalar, Varona<Scalar>>
    {
    public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Varona solver type enumeration.
       */
      using Method = enum class Method : Integer {ORDER_4 = 41, ORDER_8 = 8, ORDER_16 = 16, ORDER_32 = 32};

    private:
      Method m_method{Method::ORDER_4}; /**< Varona solver type. */

    public:

      /**
       * Class constructor for the Varona's solvers.
       */
      Varona() {}

      /**
       * Get the Varona solver name.
       * \return The Varona solver name.
       */
      constexpr std::string name_impl() const
      {
        std::ostringstream os;
        os << "Varona";
        if (this->m_method == Method::ORDER_4) {os << "4";}
        else if (this->m_method == Method::ORDER_8) {os << "8";}
        else if (this->m_method == Method::ORDER_16) {os << "16";}
        else if (this->m_method == Method::ORDER_32) {os << "32";}
        else {os << "UNKNOWN";}
        return os.str();
      }

      /**
       * Get the enumeration type of the Varona solver method.
       * \return The Varona solver enumeration type.
       */
      Method method() const {return this->m_method;}

      /**
       * Set the enumeration type of the Varona solver method.
       * \param[in] t_method The Varona solver enumeration type.
       */
      void method(Method t_method) {this->m_method = t_method;}

      /**
       * Enable the Varona's 4-th order method.
       */
      void enable_4th_order_method() {this->m_method = Method::ORDER_4;}

      /**
       * Enable the Varona's 8-th order method.
       */
      void enable_8th_order_method() {this->m_method = Method::ORDER_8;}

      /**
       * Enable the Varona's 16-th order method.
       */
      void enable_16th_order_method() {this->m_method = Method::ORDER_16;}

      /**
       * Enable the Varona's 32-th order method.
       */
      void enable_32th_order_method() {this->m_method = Method::ORDER_32;}

      /**
       * Set the Varona solver type.
       * \param[in] t_method The Varona solver type enumeration.
       */
      void set_method(Method t_method) {this->m_method = t_method;}

      /**
       * Solve the nonlinear equation \f$ f(x) = 0 \f$, with \f$ f: \mathbb{R} \rightarrow \mathbb{R} \f$.
       * \tparam FunctionLambda Function lambda type.
       * \tparam FirstDerivativeLambda First derivative lambda type.
       * \param[in] function Function lambda.
       * \param[in] first_derivative First derivative lambda.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename FirstDerivativeLambda>
      bool solve_impl(FunctionLambda && function, FirstDerivativeLambda && first_derivative, Scalar x_ini,
        Scalar & x_sol)
      {
        #define CMD "Optimist::RootFinder::Varona::solve(...): "

        // Reset internal variables
        this->reset_counters();

        // Print header
        if (this->m_verbose) {this->header();}

        // Initialize variables
        bool damped, success;
        Scalar residuals, step_norm;
        Scalar x_old, x_new, function_old, function_new, step_old, step_new;
        Scalar first_derivative_old;

        // Set initial iteration
        x_old = x_ini;
        success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old, function_old);
        OPTIMIST_ASSERT(success,
          CMD "function evaluation failed at the initial point.");

        // Algorithm iterations
        Scalar tolerance_residuals{this->m_tolerance};
        Scalar tolerance_step_norm{this->m_tolerance * this->m_tolerance};
        for (this->m_iterations = 1; this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Evaluate first derivative
          success = this->evaluate_first_derivative(std::forward<FirstDerivativeLambda>(first_derivative), x_old, first_derivative_old);
          OPTIMIST_ASSERT(success,
            CMD "first derivative evaluation failed at iteration " << this->m_iterations << ".");

          // Calculate step
          if (std::abs(first_derivative_old) < EPSILON_LOW) {
            OPTIMIST_WARNING( CMD "singular first derivative detected.");
            first_derivative_old = (first_derivative_old > 0) ? EPSILON_LOW : -EPSILON_LOW;
          }

          this->compute_step(std::forward<FunctionLambda>(function), x_old, function_old, first_derivative_old, step_old);

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
            damped = this->damp(std::forward<FunctionLambda>(function), x_old, function_old, step_old, x_new, function_new, step_new);
            OPTIMIST_ASSERT_WARNING(damped, CMD "damping failed.");
          } else {
            // Update point
            x_new = x_old + step_old;
            success = this->evaluate_function(std::forward<FunctionLambda>(function), x_new, function_new);
            OPTIMIST_ASSERT(success,
              CMD "function evaluation failed at iteration " << this->m_iterations << ".");
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

    protected:

      /**
       * Compute the step using the Varona's methods.
       * \param[in] function Function lambda.
       * \param[in] x_old Old point.
       * \param[in] function_old Old function value.
       * \param[in] first_derivative_old Old first derivative value.
       * \param[out] step_old Old step.
       */
      template <typename FunctionLambda>
      void compute_step(FunctionLambda && function, Scalar x_old, Scalar function_old, Scalar first_derivative_old,
        Scalar & step_old)
      {
        #define CMD "Optimist::RootFinder::Varona::compute_step(...): "

        bool success;
        Scalar function_y, function_z, function_w, function_h, step_tmp, t, s, u, v;
        Scalar tolerance_residuals{this->m_tolerance};

        // Base step
        step_old = -function_old/first_derivative_old;

        // Order 4 step
        if (this->m_method == Method::ORDER_4 || this->m_method == Method::ORDER_8 ||
            this->m_method == Method::ORDER_16 || this->m_method == Method::ORDER_32) {
          success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old+step_old, function_y);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed during order 4 step.");
          if (std::abs(function_y) < tolerance_residuals) {return;}
          t = function_y/function_old;
          step_tmp = this->Q(t) * (function_y/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        // Order 8 step (continued order 4)
        if (this->m_method == Method::ORDER_8 || this->m_method == Method::ORDER_16 ||
            this->m_method == Method::ORDER_32) {
          success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old+step_old, function_z);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed during order 8 step.");
          if (std::abs(function_z) < tolerance_residuals) {return;}
          s = function_z/function_y;
          step_tmp = this->W(t, s) * (function_z/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        // Order 16 step (continued order 8)
        if (this->m_method == Method::ORDER_16 || this->m_method == Method::ORDER_32) {
          success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old+step_old, function_w);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed during order 16 step.");
          if (std::abs(function_w) < tolerance_residuals) {return;}
          u = function_w/function_z;
          step_tmp = this->H(t, s, u) * (function_w/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        // Order 32 step (continued order 16)
        if (this->m_method == Method::ORDER_32) {
          success = this->evaluate_function(std::forward<FunctionLambda>(function), x_old+step_old, function_h);
          OPTIMIST_ASSERT(success,
            CMD "function evaluation failed during order 32 step.");
          if (std::abs(function_h) < tolerance_residuals) {return;}
          v = (function_h/function_w);
          step_tmp = this->J(t, s, u, v) * (function_h/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        OPTIMIST_ASSERT(std::isfinite(step_old),
          CMD "step is not finite.");

        #undef CMD
      }

      /**
       * Compute the \f$ Q \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \return The \f$ Q \f$ function value.
      **/
      static Scalar Q(Scalar t) {return 1.0 + 2.0*t;}

      /**
       * Compute the \f$ W \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \param[in] s Input value.
       * \return The \f$ W \f$ function value.
      **/
      static Scalar W(Scalar t, Scalar s) {return t*t*(1.0 - 4.0*t) + (4.0*s + 2.0)*t + s + 1.0;}

      /**
       * Compute the \f$ H \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \param[in] s Input value.
       * \param[in] u Input value.
       * \return The \f$ H \f$ function value.
      **/
      static Scalar H(Scalar t, Scalar s, Scalar u) {
        Scalar t1{t*t};
        Scalar t2{t1*t1};
        Scalar t8{s*s};
        Scalar t17{s*t8};
        Scalar t23{2.0*u};
        return
          ((8.0*u + 6.0*t2 + 4.0)*s -
          (6.0*t8 + 4.0*(s + u + 1.0))*t1 +
          2.0*t8 - 4.0*t17 + t23 + 2.0)*t +
          t1*(t8 + s + u + 1.0) +
          (1.0 - 3.0*t2 + t23)*s +
          u - t17 + 1.0;
      }

      /**
       * Compute the \f$ J \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \param[in] s Input value.
       * \param[in] u Input value.
       * \param[in] v Input value.
       * \return The \f$ J \f$ function value.
      **/
      static Scalar J(Scalar t, Scalar s, Scalar u, Scalar v) {
        Scalar t1{s*s};
        Scalar t2{t1*t1};
        Scalar t17{t*t};
        Scalar t22{u*u};
        Scalar t32{t17*t17};
        Scalar t34{t*t32};
        Scalar t37{t*t17};
        Scalar t46{1.0 + v};
        Scalar t65{u + 1.0 + v};
        Scalar t76{(-2.0*t22 + u + 4.0*v + 2.0)*u};
        return
          (-1.0 + 2.0*t)*(2.0 + 5.0*t)*u*t*t2 +
          (4.0*t + 1.0)*u*s*t2 +
          (u*t22 - 2.0*u*v - u - v - 1.0)*(4.0*t17 + 3.0*t + 1.0)*( - 1.0 + t) -
          8.0*(t22*(t17/2.0 - 1.0/4.0) + u*(t17*t32 - 5.0/8.0*t34 - 3.0/4.0*t32 +
          3.0/8.0*t37 + 3.0/4.0*t17 - t/8.0 - 1.0/4.0) + 3.0/4.0*t46*(t + 1.0/2.0)*(t - 2.0/3.0))*t*t1 +
          4.0*(t22*( - 3.0/2.0*t - 1.0/4.0) + u*(t34 - t32 - 3.0/2.0*t37 + t17/4.0 - t - 1.0/4.0) -
          t46*(t + 1.0/4.0))*s*t1 + (1.0 + v + t65*t17 - 4.0*t65*t37 - 3.0*t65*t32 + 6.0*t65*t34 +
          t76 + 4.0*(1.0 + v + t76)*t)*s;
      }

    }; // class Varona

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_VARONA4_HH
