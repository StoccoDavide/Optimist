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

#ifndef OPTIMIST_SCALAR_ROOT_FINDER_VARONA4_HXX
#define OPTIMIST_SCALAR_ROOT_FINDER_VARONA4_HXX

namespace Optimist
{
  namespace ScalarRootFinder
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
    * \includedoc docs/markdown/ScalarRootFinder/Varona.md
    */
    class Varona : public ScalarRootFinder<Varona>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{true};
      static constexpr bool requires_second_derivative{false};

      // Function types
      using Method = enum class Method : Integer {ORDER_4 = 41, ORDER_8 = 8, ORDER_16 = 16, ORDER_32 = 32}; /**< Varona solver type. */
      using FunctionWrapper         = typename ScalarRootFinder<Varona>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename ScalarRootFinder<Varona>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename ScalarRootFinder<Varona>::SecondDerivativeWrapper;

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
      std::string name_impl() const
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
      * \param[in] function Function wrapper.
      * \param[in] first_derivative First derivative wrapper.
      * \param[in] x_ini Initialization point.
      * \param[out] x_sol Solution point.
      * \return The convergence boolean flag.
      */
      bool solve_impl(FunctionWrapper function, FirstDerivativeWrapper first_derivative, Real x_ini,
        Real &x_sol)
      {
        #define CMD "Optimist::ScalarRootFinder::Varona::solve(...): "

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
        for (this->m_iterations = Integer(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
        {
          // Store trace
          this->store_trace(x_old);

          // Evaluate first derivative
          this->evaluate_first_derivative(first_derivative, x_old, first_derivative_old);

          // Calculate step
          if (std::abs(first_derivative_old) < EPSILON_LOW) {
            OPTIMIST_WARNING( CMD "singular first derivative detected.");
            first_derivative_old = (first_derivative_old > Real(0.0)) ? EPSILON_LOW : -EPSILON_LOW;
          }

          this->compute_step(function, x_old, function_old, first_derivative_old, step_old);

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

    protected:

      /**
      * Compute the step using the Varona's methods.
      * \param[in] function Function wrapper.
      * \param[in] x_old Old point.
      * \param[in] function_old Old function value.
      * \param[in] first_derivative_old Old first derivative value.
      * \param[out] step_old Old step.
      */
      void compute_step(FunctionWrapper function, Real x_old, Real function_old, Real first_derivative_old,
        Real & step_old)
      {
        Real function_y, function_z, function_w, function_h, step_tmp, t, s, u, v;
        Real tolerance_residuals{this->m_tolerance};
        //Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};

        // Base step
        step_old = -function_old/first_derivative_old;

        // Order 4 step
        if (this->m_method == Method::ORDER_4 || this->m_method == Method::ORDER_8 ||
            this->m_method == Method::ORDER_16 || this->m_method == Method::ORDER_32) {
          this->evaluate_function(function, x_old+step_old, function_y);
          if (std::abs(function_y) < tolerance_residuals) {return;}
          t = function_y/function_old;
          step_tmp = this->Q(t) * (function_y/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        // Order 8 step (continued order 4)
        if (this->m_method == Method::ORDER_8 || this->m_method == Method::ORDER_16 ||
            this->m_method == Method::ORDER_32) {
          this->evaluate_function(function, x_old+step_old, function_z);
          if (std::abs(function_z) < tolerance_residuals) {return;}
          s = function_z/function_y;
          step_tmp = this->W(t, s) * (function_z/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        // Order 16 step (continued order 8)
        if (this->m_method == Method::ORDER_16 || this->m_method == Method::ORDER_32) {
          this->evaluate_function(function, x_old+step_old, function_w);
          if (std::abs(function_w) < tolerance_residuals) {return;}
          u = function_w/function_z;
          step_tmp = this->H(t, s, u) * (function_w/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        // Order 32 step (continued order 16)
        if (this->m_method == Method::ORDER_32) {
          this->evaluate_function(function, x_old+step_old, function_h);
          if (std::abs(function_h) < tolerance_residuals) {return;}
          v = (function_h/function_w);
          step_tmp = this->J(t, s, u, v) * (function_h/first_derivative_old);
          if (std::isfinite(step_tmp)) {step_old -= step_tmp;}
          else {return;}
        }

        OPTIMIST_ASSERT(std::isfinite(step_old),
          "Optimist::ScalarRootFinder::Varona::compute_step(...): step is not finite.");
      }

      /**
      * Compute the \f$ Q \f$ function for the Varona's methods.
      * \param[in] t Input value.
      * \return The \f$ Q \f$ function value.
      **/
      static Real Q(Real t) {return 1.0 + 2.0*t;}

      /**
      * Compute the \f$ W \f$ function for the Varona's methods.
      * \param[in] t Input value.
      * \param[in] s Input value.
      * \return The \f$ W \f$ function value.
      **/
      static Real W(Real t, Real s) {return t*t*(1.0 - 4.0*t) + (4.0*s + 2.0)*t + s + 1.0;}

      /**
      * Compute the \f$ H \f$ function for the Varona's methods.
      * \param[in] t Input value.
      * \param[in] s Input value.
      * \param[in] u Input value.
      * \return The \f$ H \f$ function value.
      **/
      static Real H(Real t, Real s, Real u) {
        Real t1{t*t};
        Real t2{t1*t1};
        Real t8{s*s};
        Real t17{s*t8};
        Real t23{2.0*u};
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
      static Real J(Real t, Real s, Real u, Real v) {
        Real t1{s*s};
        Real t2{t1*t1};
        Real t17{t*t};
        Real t22{u*u};
        Real t32{t17*t17};
        Real t34{t*t32};
        Real t37{t*t17};
        Real t46{1.0 + v};
        Real t65{u + 1.0 + v};
        Real t76{(-2.0*t22 + u + 4.0*v + 2.0)*u};
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

  } // namespace ScalarRootFinder

} // namespace Optimist

#endif // OPTIMIST_SCALAR_ROOT_FINDER_VARONA4_HXX
