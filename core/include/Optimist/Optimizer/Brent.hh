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

#ifndef OPTIMIST_OPTIMIZER_BRENT_HH
#define OPTIMIST_OPTIMIZER_BRENT_HH

#include "Optimist/Optimizer/ScalarMinimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for Brent's minimization method.
     *
     * \includedoc docs/markdown/Optimizer/Brent.md
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class Brent : public ScalarMinimizer<Scalar, Brent<Scalar>> {
     public:
      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Class constructor for the Brent solver.
       */
      Brent() {}

      /**
       * Get the Brent solver name.
       */
      constexpr std::string name_impl() const {
        return "Brent";
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

        // x: best point, w: second best, v: previous value of w
        const Scalar cgold{1 - 1 / this->GOLDEN};
        Scalar x{b}, w{b}, v{b}, fx{fb}, fw{fb}, fv{fb};
        Scalar d{0}, e{0};
        std::string notes{"Initialize"};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          const Scalar xm{(a + c) / 2};
          const Scalar tol1{this->m_tolerance * (1 + std::abs(x))};
          const Scalar tol2{2 * tol1};
          if (this->m_verbose) {
            this->info(c - a, notes);
          }
          if (std::abs(x - xm) <= tol2 - (c - a) / 2) {
            this->m_converged = true;
            break;
          }

          bool golden{true};
          if (std::abs(e) > tol1) {
            // Trial parabolic fit through x, v, w
            const Scalar r{(x - w) * (fx - fv)};
            Scalar q{(x - v) * (fx - fw)};
            Scalar p{(x - v) * q - (x - w) * r};
            q = 2 * (q - r);
            if (q > 0) {
              p = -p;
            }
            q = std::abs(q);
            const Scalar e_old{e};
            e = d;
            if (std::abs(p) < std::abs(q * e_old / 2) && p > q * (a - x) &&
                p < q * (c - x)) {
              d = p / q;
              const Scalar u{x + d};
              if (u - a < tol2 || c - u < tol2) {
                d = xm >= x ? tol1 : -tol1;
              }
              golden = false;
              notes  = "Parabolic";
            }
          }
          if (golden) {
            e     = x >= xm ? a - x : c - x;
            d     = cgold * e;
            notes = "Golden";
          }

          const Scalar u{std::abs(d) >= tol1 ? x + d
                                             : x + (d >= 0 ? tol1 : -tol1)};
          const Scalar fu{this->eval(function, u)};
          if (fu <= fx) {
            if (u >= x) {
              a = x;
            } else {
              c = x;
            }
            v  = w;
            fv = fw;
            w  = x;
            fw = fx;
            x  = u;
            fx = fu;
          } else {
            if (u < x) {
              a = u;
            } else {
              c = u;
            }
            if (fu <= fw || w == x) {
              v  = w;
              fv = fw;
              w  = u;
              fw = fu;
            } else if (fu <= fv || v == x || v == w) {
              v  = u;
              fv = fu;
            }
          }
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = x;
        return this->m_converged;
      }

    };  // class Brent

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_BRENT_HH
