/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Optimist project is distributed under the BSD 2-Clause License.       *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef OPTIMIST_OPTIMIZER_SCALAR_MINIMIZER_HH
#define OPTIMIST_OPTIMIZER_SCALAR_MINIMIZER_HH

#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Base class for derivative-free minimizers of functions of one
     * variable.
     *
     * The search starts from a bracketing triplet \f$ a < b < c \f$ (or \f$ a
     * > b > c \f$) with \f$ f(b) < f(a) \f$ and \f$ f(b) \leq f(c) \f$. If
     * finite bounds are set, the bracket is the whole bounded interval (the
     * function is assumed unimodal there); otherwise it is found by golden
     * expansion downhill from the initial guess.
     *
     * \tparam Scalar Floating-point number type.
     * \tparam DerivedSolver Derived solver class.
     */
    template <typename Scalar, typename DerivedSolver>
      requires TypeTrait<Scalar>::IsScalar
    class ScalarMinimizer : public Optimizer<Scalar, DerivedSolver> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Golden ratio \f$ \varphi = (1 + \sqrt{5})/2 \f$.
       */
      inline static const Scalar GOLDEN{(1 + std::sqrt(Scalar(5))) / 2};

     protected:
      Scalar m_step{1.0};           /**< Initial bracketing step. */
      Integer m_max_expansions{60}; /**< Maximum bracketing expansions. */

     public:
      /**
       * Get the initial bracketing step.
       */
      Scalar step() const {
        return this->m_step;
      }

      /**
       * Set the initial bracketing step.
       */
      void step(const Scalar t_step) {
        OPTIMIST_ASSERT(t_step > 0,
                        "Optimist::Optimizer::ScalarMinimizer::step(...): "
                        "invalid input detected.");
        this->m_step = t_step;
      }

     protected:
      /**
       * Evaluate the function, mapping failures and non-finite values to
       * \f$ +\infty \f$.
       */
      template <typename FunctionLambda>
      Scalar eval(FunctionLambda &function, const Scalar x) {
        Scalar f;
        if (!this->evaluate_function(function, x, f) || !std::isfinite(f)) {
          return INFTY;
        }
        return f;
      }

      /**
       * Find a bracketing triplet.
       * \return True if a bracket was found.
       */
      template <typename FunctionLambda>
      bool bracket(FunctionLambda &function,
                   const Scalar x_ini,
                   Scalar &a,
                   Scalar &b,
                   Scalar &c,
                   Scalar &fa,
                   Scalar &fb,
                   Scalar &fc) {
        if (this->has_finite_bounds()) {
          a  = this->m_lower_bound;
          c  = this->m_upper_bound;
          b  = std::clamp(x_ini, a, c);
          if (b == a || b == c) {
            b = a + (c - a) / (1 + GOLDEN);
          }
          fa = this->eval(function, a);
          fb = this->eval(function, b);
          fc = this->eval(function, c);
          return true;
        }
        a  = x_ini;
        b  = x_ini + this->m_step;
        fa = this->eval(function, a);
        fb = this->eval(function, b);
        if (fb > fa) {
          std::swap(a, b);
          std::swap(fa, fb);
        }
        c  = b + GOLDEN * (b - a);
        fc = this->eval(function, c);
        for (Integer i{0}; i < this->m_max_expansions && !(fb <= fc); ++i) {
          a  = b;
          fa = fb;
          b  = c;
          fb = fc;
          c  = b + GOLDEN * (b - a);
          fc = this->eval(function, c);
          if (!std::isfinite(c)) {
            return false;
          }
        }
        return fb <= fc && fb <= fa;
      }

    };  // class ScalarMinimizer

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_SCALAR_MINIMIZER_HH
