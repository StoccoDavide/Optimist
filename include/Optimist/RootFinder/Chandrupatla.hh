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

#ifndef OPTIMIST_ROOTFINDER_CHANDRUPATLA_HH
#define OPTIMIST_ROOTFINDER_CHANDRUPATLA_HH

#include "Optimist/RootFinder/Bracketing.hh"

namespace Optimist
{
  namespace RootFinder
  {

    /*\
     |    ____ _                     _                        _   _
     |   / ___| |__   __ _ _ __   __| |_ __ _   _ _ __   __ _| |_| | __ _
     |  | |   | '_ \ / _` | '_ \ / _` | '__| | | | '_ \ / _` | __| |/ _` |
     |  | |___| | | | (_| | | | | (_| | |  | |_| | |_) | (_| | |_| | (_| |
     |   \____|_| |_|\__,_|_| |_|\__,_|_|   \__,_| .__/ \__,_|\__|_|\__,_|
     |                                           |_|
    \*/

    /**
    * \brief Class container for the Chandrupatla algorithm.
    *
    * The Chandrupatla algorithm allows to find the roots of a scalar function \f$f(x)\f$ in a given
    * interval \f$[a, b]\f$. The algorithm is based on the work of T. Chandrupatla, *A new hybrid
    * quadratic/bisection algorithm for finding the zero of a nonlinear function without using
    * derivatives*, Advances in Engineering Software, 28 (1997), pp. 145-149.
    * \tparam Real Scalar number type.
    */
    template <typename Real>
    class Chandrupatla : public Bracketing<Real, Chandrupatla<Real>>
    {
    public:
      static constexpr bool requires_function{true};
      static constexpr bool requires_first_derivative{false};
      static constexpr bool requires_second_derivative{false};

      OPTIMIST_BASIC_CONSTANTS(Real) /**< Basic constants. */

      // Function types
      using FunctionWrapper         = typename RootFinder<Real, 1, Chandrupatla>::FunctionWrapper;
      using FirstDerivativeWrapper  = typename RootFinder<Real, 1, Chandrupatla>::FirstDerivativeWrapper;
      using SecondDerivativeWrapper = typename RootFinder<Real, 1, Chandrupatla>::SecondDerivativeWrapper;

      /**
      * Class constructor for the Algorithm 748.
      */
      Chandrupatla() {}

      /**
      * Get the Algorithm 748 solver name.
      * \return The Algorithm 748 solver name.
      */
      std::string name_impl() const {return "Chandrupatla";}

      /**
      * Finds either an exact solution or an approximate solution of the equation \f$f(x) = 0\f$ in
      * the interval \f$[a, b]\f$.
      * \param[in] function Function wrapper.
      * \return The approximate root.
      */
      Real find_root_impl(FunctionWrapper function) {

        #define CMD "Optimist::ScalarRootfinder::Chandrupatla::find_root_impl(...): "

        Real tolerance_step{this->m_tolerance_bracketing};
        Real tolerance_function{this->m_tolerance_bracketing};

        // Check for trivial solution
        this->m_converged = this->m_fa == 0; if (this->m_converged) {return this->m_a;}
        this->m_converged = this->m_fb == 0; if (this->m_converged) {return this->m_b;}

        // Initialize to dumb values
        this->m_e  = QUIET_NAN;
        this->m_fe = QUIET_NAN;

        // While f(left) or f(right) are infinite perform bisection
        while (!(std::isfinite(this->m_fa) && std::isfinite(this->m_fb))) {
          ++this->m_iterations;
          this->m_c  = (this->m_a + this->m_b)/2.0;
          this->evaluate_function(function, this->m_c, this->m_fc);
          this->m_converged = this->m_fc == 0;
          if (this->m_converged) {return this->m_c;}
          if (this->m_fa*this->m_fc < 0.0) {
            // -> [a, c]
            this->m_b = this->m_c; this->m_fb = this->m_fc;
          } else {
            // -> [c, b]
            this->m_a = this->m_c; this->m_fa = this->m_fc;
          }

          // Check for convergence
          Real abs_fa{std::abs(this->m_fa)};
          Real abs_fb{std::abs(this->m_fb)};
          this->m_converged = (this->m_b - this->m_a) < tolerance_step ||
                               abs_fa < tolerance_function ||abs_fb < tolerance_function;
          if (this->m_converged) {return abs_fb < abs_fa ? this->m_b : this->m_a;}
        }

        Real t{0.5};

        while (++this->m_iterations < this->m_max_iterations) {

          Real direction{this->m_b - this->m_a};
          Real tolerance{tolerance_step/(2.0*std::abs(direction))};
          this->m_converged = tolerance >= 0.5;
          if (this->m_converged) {break;}

          Real abs_fa{std::abs(this->m_fa)};
          Real abs_fb{std::abs(this->m_fb)};
          this->m_converged = abs_fa < tolerance_function || abs_fb < tolerance_function;
          if (this->m_converged) {break;}

          if (t < tolerance) {t = tolerance;}
          else if (t > 1.0 - tolerance) {t = 1.0 - tolerance;}

          Real c{this->m_a + t * direction};
          Real fc; this->evaluate_function(function, c, fc);
          if (this->m_verbose) {this->info(fc);}

          OPTIMIST_ASSERT(std::isfinite(fc), CMD "function is not finite.");

          this->m_converged = fc == 0;
          if (this->m_converged) {
            this->m_a = this->m_b = c; this->m_fa = this->m_fb = 0;
            break;
          }

          // Arrange 2-1-3: 2-1 interval, 1 middle, 3 discarded point
          Real d, fd;
          if ((0.0 < fc) == (0.0 < this->m_fa)) {
            // a -> d  --> [d,a,b] = [a,c,b]
            d = this->m_a; fd = this->m_fa;
          } else {
            // b -> d, a -> b  --> [b,a,d] = [a,c,b]
            d = this->m_b; fd = this->m_fb;
            this->m_b = this->m_a; this->m_fb = this->m_fa;
          }
          //  c => a
          this->m_a = c; this->m_fa = fc;

          // If inverse quadratic interpolation holds use it
          Real ba{this->m_b - this->m_a};
          Real fba{this->m_fb - this->m_fa};
          Real bd{this->m_b - d};
          Real fbd{this->m_fb - fd};

          Real xi{ba/bd};
          Real ph{fba/fbd};
          Real fl{1.0 - std::sqrt(1.0 - xi)};
          Real fh{std::sqrt(xi)};

          if (fl < ph && ph < fh) {
            Real da{d - this->m_a};
            Real fda{fd - this->m_fa};
            t = (this->m_fa/fba) * (fd/fbd) - (this->m_fa/fda) * (this->m_fb/fbd) * (da/ba);
          } else {
            t = 0.5;
          }
        }
        if (this->m_a > this->m_b) {std::swap(this->m_a, this->m_b); std::swap(this->m_fa, this->m_fb);}
        return std::abs(this->m_fa) < std::abs(this->m_fb) ? this->m_a : this->m_b;;

        #undef CMD
      }

    }; // class Chandrupatla

  } // namespace RootFinder

} // namespace Optimist

#endif // OPTIMIST_ROOTFINDER_CHANDRUPATLA_HH
