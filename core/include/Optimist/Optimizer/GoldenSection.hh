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

#ifndef OPTIMIST_OPTIMIZER_GOLDEN_SECTION_HH
#define OPTIMIST_OPTIMIZER_GOLDEN_SECTION_HH

#include "Optimist/Optimizer/ScalarMinimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the golden-section search.
     *
     * \includedoc docs/markdown/Optimizer/GoldenSection.md
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class GoldenSection
        : public ScalarMinimizer<Scalar, GoldenSection<Scalar>> {
     public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the golden-section solver.
       */
      GoldenSection() {}

      /**
       * Get the golden-section solver name.
       */
      constexpr std::string name_impl() const {
        return "GoldenSection";
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda>
      bool solve_impl(FunctionLambda &&function,
                      const Scalar x_ini,
                      Scalar &x_sol) {
        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        Scalar a, b, c, fa, fb, fc;
        if (!this->bracket(function, x_ini, a, b, c, fa, fb, fc)) {
          if (this->m_verbose) {
            this->bottom();
          }
          x_sol = fa < fc ? a : c;
          return false;
        }
        if (a > c) {
          std::swap(a, c);
          std::swap(fa, fc);
        }

        // Two interior points in golden ratio
        const Scalar r{1 / this->GOLDEN};
        Scalar x1{c - r * (c - a)}, x2{a + r * (c - a)};
        Scalar f1{this->eval(function, x1)}, f2{this->eval(function, x2)};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          const Scalar width{c - a};
          if (this->m_verbose) {
            this->info(width);
          }
          if (width <= this->m_tolerance * (1 + std::abs(x1) + std::abs(x2))) {
            this->m_converged = true;
            break;
          }
          if (f1 < f2) {
            c  = x2;
            x2 = x1;
            f2 = f1;
            x1 = c - r * (c - a);
            f1 = this->eval(function, x1);
          } else {
            a  = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + r * (c - a);
            f2 = this->eval(function, x2);
          }
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = f1 < f2 ? x1 : x2;
        return this->m_converged;
      }

    };  // class GoldenSection

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_GOLDEN_SECTION_HH
