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

#ifndef OPTIMIST_ROOTFINDER_CHANDRUPATLA_HH
#define OPTIMIST_ROOTFINDER_CHANDRUPATLA_HH

#include "Optimist/RootFinder/Bracketing.hh"

namespace Optimist {
  namespace RootFinder {

    /**
     * \brief Class container for the Chandrupatla algorithm.
     *
     * The Chandrupatla algorithm allows to find the roots of a scalar function
     * \f$f(x)\f$ in a given interval \f$[a, b]\f$. The algorithm is based on
     * the work of T. Chandrupatla, *A new hybrid quadratic/bisection algorithm
     * for finding the zero of a nonlinear function without using derivatives*,
     * Advances in Engineering Software, 28 (1997), pp. 145-149.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    class Chandrupatla : public Bracketing<Scalar, Chandrupatla<Scalar>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Chandrupatla solver.
       */
      Chandrupatla() {}

      /**
       * Get the Chandrupatla solver name.
       * \return The Chandrupatla solver name.
       */
      constexpr std::string name_impl() const {
        return "Chandrupatla";
      }

      /**
       * Finds either an exact solution or an approximate solution of the
       * equation
       * \f$f(x) = 0\f$ in the interval \f$[a, b]\f$.
       * \tparam FunctionLambda Function lambda type.
       * \param[in] function Function lambda.
       * \return The approximate root.
       */
      template <typename FunctionLambda>
      Scalar find_root_impl(FunctionLambda &&function) {
#define CMD "Optimist::RootFinder::Chandrupatla::find_root_impl(...): "

        bool success;
        Scalar tolerance_step{this->m_tolerance_bracketing};
        Scalar tolerance_function{this->m_tolerance_bracketing};

        // Check for trivial solution
        this->m_converged = this->m_fa == 0;
        if (this->m_verbose) {
          this->info(std::abs(this->m_fa));
        }
        if (this->m_converged) {
          return this->m_a;
        }
        this->m_converged = this->m_fb == 0;
        if (this->m_verbose) {
          this->info(std::abs(this->m_fb));
        }
        if (this->m_converged) {
          return this->m_b;
        }

        // Initialize to dumb values
        this->m_e  = QUIET_NAN;
        this->m_fe = QUIET_NAN;

        // While f(left) or f(right) are infinite perform bisection
        while (!(std::isfinite(this->m_fa) && std::isfinite(this->m_fb))) {
          ++this->m_iterations;
          this->m_c = (this->m_a + this->m_b) / 2.0;
          success =
              this->evaluate_function(std::forward<FunctionLambda>(function),
                                      this->m_c,
                                      this->m_fc);
          OPTIMIST_ASSERT(success,
                          CMD "function evaluation failed at iteration "
                              << this->m_iterations << ".");
          this->m_converged = this->m_fc == 0;
          if (this->m_verbose) {
            this->info(std::abs(this->m_fc));
          }
          if (this->m_converged) {
            return this->m_c;
          }
          if (this->m_fa * this->m_fc < 0) {
            // -> [a, c]
            this->m_b  = this->m_c;
            this->m_fb = this->m_fc;
          } else {
            // -> [c, b]
            this->m_a  = this->m_c;
            this->m_fa = this->m_fc;
          }

          // Check for convergence
          Scalar abs_fa{std::abs(this->m_fa)};
          Scalar abs_fb{std::abs(this->m_fb)};
          if (this->m_verbose) {
            this->info(abs_fa < abs_fb ? abs_fa : abs_fb);
          }
          this->m_converged = (this->m_b - this->m_a) < tolerance_step ||
                              abs_fa < tolerance_function ||
                              abs_fb < tolerance_function;
          if (this->m_converged) {
            return abs_fb < abs_fa ? this->m_b : this->m_a;
          }
        }

        Scalar t{0.5};

        while (++this->m_iterations < this->m_max_iterations) {
          Scalar direction{this->m_b - this->m_a};
          Scalar tolerance{tolerance_step /
                           (static_cast<Scalar>(2.0) * std::abs(direction))};
          this->m_converged = tolerance >= 0.5;
          if (this->m_converged) {
            break;
          }

          Scalar abs_fa{std::abs(this->m_fa)};
          Scalar abs_fb{std::abs(this->m_fb)};
          if (this->m_verbose) {
            this->info(abs_fa < abs_fb ? abs_fa : abs_fb);
          }
          this->m_converged =
              abs_fa < tolerance_function || abs_fb < tolerance_function;
          if (this->m_converged) {
            break;
          }

          if (t < tolerance) {
            t = tolerance;
          } else if (t > 1.0 - tolerance) {
            t = 1.0 - tolerance;
          }

          Scalar c{this->m_a + t * direction};
          Scalar fc;
          success =
              this->evaluate_function(std::forward<FunctionLambda>(function),
                                      c,
                                      fc);
          OPTIMIST_ASSERT(success,
                          CMD "function evaluation failed at iteration "
                              << this->m_iterations << ".");
          if (this->m_verbose) {
            this->info(fc);
          }

          OPTIMIST_ASSERT(std::isfinite(fc), CMD "function is not finite.");

          this->m_converged = fc == 0;
          if (this->m_converged) {
            this->m_a = this->m_b = c;
            this->m_fa = this->m_fb = 0;
            break;
          }

          // Arrange 2-1-3: 2-1 interval, 1 middle, 3 discarded point
          Scalar d, fd;
          if ((0.0 < fc) == (0.0 < this->m_fa)) {
            // a -> d  --> [d,a,b] = [a,c,b]
            d  = this->m_a;
            fd = this->m_fa;
          } else {
            // b -> d, a -> b  --> [b,a,d] = [a,c,b]
            d          = this->m_b;
            fd         = this->m_fb;
            this->m_b  = this->m_a;
            this->m_fb = this->m_fa;
          }
          //  c => a
          this->m_a  = c;
          this->m_fa = fc;

          // If inverse quadratic interpolation holds use it
          Scalar ba{this->m_b - this->m_a};
          Scalar fba{this->m_fb - this->m_fa};
          Scalar bd{this->m_b - d};
          Scalar fbd{this->m_fb - fd};

          Scalar xi{ba / bd};
          Scalar ph{fba / fbd};
          Scalar fl{static_cast<Scalar>(1.0) -
                    std::sqrt(static_cast<Scalar>(1.0) - xi)};
          Scalar fh{std::sqrt(xi)};

          if (fl < ph && ph < fh) {
            Scalar da{d - this->m_a};
            Scalar fda{fd - this->m_fa};
            t = (this->m_fa / fba) * (fd / fbd) -
                (this->m_fa / fda) * (this->m_fb / fbd) * (da / ba);
          } else {
            t = 0.5;
          }
        }
        if (this->m_a > this->m_b) {
          std::swap(this->m_a, this->m_b);
          std::swap(this->m_fa, this->m_fb);
        }
        return std::abs(this->m_fa) < std::abs(this->m_fb) ? this->m_a
                                                           : this->m_b;
        ;

#undef CMD
      }

    };  // class Chandrupatla

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_CHANDRUPATLA_HH
