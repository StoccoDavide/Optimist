/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_ROOTFINDER_VARONA_HH
#define OPTIMIST_ROOTFINDER_VARONA_HH

#include "Optimist/RootFinder/OpenMethod.hh"

namespace Optimist {

  namespace RootFinder {

    /**
     * \brief Class container for Varona's family of high-order methods.
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class Varona : public OpenMethod<Scalar, Varona<Scalar>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Order of the method.
       */
      using Method = enum class Method : Integer {
        ORDER_4  = 4,  /**< Fourth order (2 function evaluations). */
        ORDER_8  = 8,  /**< Eighth order (3 function evaluations). */
        ORDER_16 = 16, /**< Sixteenth order (4 function evaluations). */
        ORDER_32 = 32  /**< Thirty-second order (5 function evaluations). */
      };

     private:
      Method m_method{Method::ORDER_4}; /**< Order of the method. */

     public:
      /**
       * Class constructor for the Varona solver.
       */
      Varona() {}

      /**
       * Get the Varona solver name.
       * \return The Varona solver name.
       */
      constexpr std::string name_impl() const {
        return "Varona";
      }

      /**
       * Get the order of the method.
       * \return The order of the method.
       */
      Method method() const {
        return this->m_method;
      }

      /**
       * Set the order of the method.
       * \param[in] t_method The order of the method.
       */
      void method(const Method t_method) {
        this->m_method = t_method;
      }

      /**
       * Enable the fourth-order method.
       */
      void enable_4th_order_method() {
        this->m_method = Method::ORDER_4;
      }

      /**
       * Enable the eighth-order method.
       */
      void enable_8th_order_method() {
        this->m_method = Method::ORDER_8;
      }

      /**
       * Enable the sixteenth-order method.
       */
      void enable_16th_order_method() {
        this->m_method = Method::ORDER_16;
      }

      /**
       * Enable the thirty-second-order method.
       */
      void enable_32th_order_method() {
        this->m_method = Method::ORDER_32;
      }

      /**
       * Varona step: a Newton step followed by up to four corrections, each
       * needing one more function evaluation and doubling the order. The
       * weight functions are expansions around \f$ t = f(y)/f(x) = 0 \f$, so
       * far from the root (\f$ |t| > 1/2 \f$) the Newton step is taken.
       * \param[in] function Function.
       * \param[in] x Current point.
       * \param[in] f Function value.
       * \param[in] df First derivative value.
       * \param[in] ddf Second derivative value (unused).
       * \return The step.
       */
      template <typename FunctionLambda>
      Scalar step_impl(FunctionLambda &function,
                       const Scalar x,
                       const Scalar f,
                       const Scalar df,
                       [[maybe_unused]] const Scalar ddf) {
#define CMD "Optimist::RootFinder::Varona::step_impl(...): "

        const Integer order{static_cast<Integer>(this->m_method)};
        const Scalar tolerance{this->m_tolerance};
        Scalar step{-f / df};

        // Evaluate the function at x + step, returning false if converged
        auto next = [&](Scalar &value) -> bool {
          OPTIMIST_ASSERT(this->evaluate_function(function, x + step, value),
                          CMD "function evaluation failed.");
          return std::abs(value) >= tolerance;
        };

        Scalar f_y, f_z, f_w, f_h;
        if (!next(f_y)) {
          return step;
        }
        const Scalar t{f_y / f};
        if (std::abs(t) > 0.5) {
          return step;
        }
        step -= Q(t) * f_y / df;
        if (order < 8 || !next(f_z)) {
          return step;
        }
        const Scalar s{f_z / f_y};
        step -= W(t, s) * f_z / df;
        if (order < 16 || !next(f_w)) {
          return step;
        }
        const Scalar u{f_w / f_z};
        step -= H(t, s, u) * f_w / df;
        if (order < 32 || !next(f_h)) {
          return step;
        }
        const Scalar w{f_h / f_w};
        return step - J(t, s, u, w) * f_h / df;

#undef CMD
      }

     private:
      /**
       * Compute the \f$ Q \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \return The \f$ Q \f$ function value.
       **/
      static Scalar Q(Scalar t) {
        return 1.0 + 2.0 * t;
      }

      /**
       * Compute the \f$ W \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \param[in] s Input value.
       * \return The \f$ W \f$ function value.
       **/
      static Scalar W(Scalar t, Scalar s) {
        return t * t * (1.0 - 4.0 * t) + (4.0 * s + 2.0) * t + s + 1.0;
      }

      /**
       * Compute the \f$ H \f$ function for the Varona's methods.
       * \param[in] t Input value.
       * \param[in] s Input value.
       * \param[in] u Input value.
       * \return The \f$ H \f$ function value.
       **/
      static Scalar H(Scalar t, Scalar s, Scalar u) {
        Scalar t1{t * t};
        Scalar t2{t1 * t1};
        Scalar t8{s * s};
        Scalar t17{s * t8};
        Scalar t23{static_cast<Scalar>(2.0) * u};
        return ((8.0 * u + 6.0 * t2 + 4.0) * s -
                (6.0 * t8 + 4.0 * (s + u + 1.0)) * t1 + 2.0 * t8 - 4.0 * t17 +
                t23 + 2.0) *
                   t +
               t1 * (t8 + s + u + 1.0) + (1.0 - 3.0 * t2 + t23) * s + u - t17 +
               1.0;
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
        Scalar t1{s * s};
        Scalar t2{t1 * t1};
        Scalar t17{t * t};
        Scalar t22{u * u};
        Scalar t32{t17 * t17};
        Scalar t34{t * t32};
        Scalar t37{t * t17};
        Scalar t46{static_cast<Scalar>(1.0) + v};
        Scalar t65{u + static_cast<Scalar>(1.0) + v};
        Scalar t76{static_cast<Scalar>(-2.0 * t22 + u + 4.0 * v + 2.0) * u};
        return (-1.0 + 2.0 * t) * (2.0 + 5.0 * t) * u * t * t2 +
               (4.0 * t + 1.0) * u * s * t2 +
               (u * t22 - 2.0 * u * v - u - v - 1.0) *
                   (4.0 * t17 + 3.0 * t + 1.0) * (-1.0 + t) -
               8.0 *
                   (t22 * (t17 / 2.0 - 1.0 / 4.0) +
                    u * (t17 * t32 - 5.0 / 8.0 * t34 - 3.0 / 4.0 * t32 +
                         3.0 / 8.0 * t37 + 3.0 / 4.0 * t17 - t / 8.0 -
                         1.0 / 4.0) +
                    3.0 / 4.0 * t46 * (t + 1.0 / 2.0) * (t - 2.0 / 3.0)) *
                   t * t1 +
               4.0 *
                   (t22 * (-3.0 / 2.0 * t - 1.0 / 4.0) +
                    u * (t34 - t32 - 3.0 / 2.0 * t37 + t17 / 4.0 - t -
                         1.0 / 4.0) -
                    t46 * (t + 1.0 / 4.0)) *
                   s * t1 +
               (1.0 + v + t65 * t17 - 4.0 * t65 * t37 - 3.0 * t65 * t32 +
                6.0 * t65 * t34 + t76 + 4.0 * (1.0 + v + t76) * t) *
                   s;
      }

    };  // class Varona

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_VARONA_HH
