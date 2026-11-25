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

#ifndef OPTIMIST_ROOTFINDER_ITP_HH
#define OPTIMIST_ROOTFINDER_ITP_HH

#include "Optimist/RootFinder/Bracketing.hh"

namespace Optimist {

  namespace RootFinder {

    /**
     * \brief Class container for the ITP (Interpolate, Truncate, Project)
     * method of Oliveira and Takahashi (2020).
     *
     * \includedoc docs/markdown/RootFinder/ITP.md
     *
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
      requires TypeTrait<Scalar>::IsScalar
    class ITP : public Bracketing<Scalar, ITP<Scalar>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{false};
      static constexpr bool RequiresSecondDerivative{false};

      OPTIMIST_BASIC_CONSTANTS(Scalar)

     private:
      Scalar m_kappa_1{0.2}; /**< Truncation factor (relative to b - a). */
      Scalar m_kappa_2{2.0}; /**< Truncation exponent in [1, 1 + phi). */
      Integer m_n_0{1};      /**< Slack iterations over bisection. */

     public:
      /**
       * Class constructor for the ITP solver.
       */
      ITP() {}

      /**
       * Get the ITP solver name.
       */
      constexpr std::string name_impl() const {
        return "ITP";
      }

      /**
       * Set the hyper-parameters of the method.
       * \param[in] t_kappa_1 Truncation factor \f$ \kappa_1 > 0 \f$, relative
       * to the initial interval width.
       * \param[in] t_kappa_2 Truncation exponent \f$ \kappa_2 \in [1, 1 +
       * \phi) \f$, with \f$ \phi \f$ the golden ratio.
       * \param[in] t_n_0 Number of iterations allowed in excess of bisection.
       */
      void parameters(const Scalar t_kappa_1,
                      const Scalar t_kappa_2,
                      const Integer t_n_0) {
        OPTIMIST_ASSERT(t_kappa_1 > 0 && t_kappa_2 >= 1 &&
                            t_kappa_2 < 1 + (1 + std::sqrt(Scalar(5))) / 2 &&
                            t_n_0 >= 0,
                        "Optimist::RootFinder::ITP::parameters(...): invalid "
                        "input detected.");
        this->m_kappa_1 = t_kappa_1;
        this->m_kappa_2 = t_kappa_2;
        this->m_n_0     = t_n_0;
      }

      /**
       * Find the root of the function in the bracket \f$ [a, b] \f$.
       * \param[in] function Function.
       * \return The root.
       */
      template <typename FunctionLambda>
      Scalar find_root_impl(FunctionLambda &&function) {
#define CMD "Optimist::RootFinder::ITP::find_root_impl(...): "

        Scalar a{this->m_a}, b{this->m_b}, fa{this->m_fa}, fb{this->m_fb};
        if (a > b) {
          std::swap(a, b);
          std::swap(fa, fb);
        }
        if (fa == 0 || fb == 0) {
          this->m_converged = true;
          return fa == 0 ? a : b;
        }

        // Orient the function so that f(a) < 0 < f(b)
        const Scalar sign{fa < 0 ? Scalar(1) : Scalar(-1)};
        fa *= sign;
        fb *= sign;

        const Scalar eps{std::max(
            this->m_tolerance_bracketing,
            EPSILON * std::max({Scalar(1), std::abs(a), std::abs(b)}))};
        const Scalar kappa_1{this->m_kappa_1 / (b - a)};
        const Integer n_half{static_cast<Integer>(
            std::ceil(std::log2((b - a) / (2 * eps))))};
        const Integer n_max{std::max(n_half, Integer(0)) + this->m_n_0};

        for (Integer j{0}; this->m_iterations < this->m_max_iterations; ++j) {
          ++this->m_iterations;
          if (this->m_verbose) {
            this->info(std::min(std::abs(fa), std::abs(fb)));
          }
          if (b - a <= 2 * eps) {
            this->m_converged = true;
            break;
          }

          // Interpolation (regula falsi)
          const Scalar x_half{(a + b) / 2};
          const Scalar x_f{(fb * a - fa * b) / (fb - fa)};

          // Truncation
          const Scalar sigma{x_half >= x_f ? Scalar(1) : Scalar(-1)};
          const Scalar delta{kappa_1 * std::pow(b - a, this->m_kappa_2)};
          const Scalar x_t{delta <= std::abs(x_half - x_f) ? x_f + sigma * delta
                                                           : x_half};

          // Projection onto the minmax interval
          const Scalar r{eps * std::ldexp(Scalar(1), n_max - j) - (b - a) / 2};
          const Scalar x_itp{std::abs(x_t - x_half) <= r ? x_t
                                                         : x_half - sigma * r};

          // Update the bracket
          Scalar f_itp;
          bool success{this->evaluate_function(function, x_itp, f_itp)};
          OPTIMIST_ASSERT(success && std::isfinite(f_itp),
                          CMD "function evaluation failed at iteration "
                              << this->m_iterations << ".");
          f_itp *= sign;
          if (f_itp > 0) {
            b  = x_itp;
            fb = f_itp;
          } else if (f_itp < 0) {
            a  = x_itp;
            fa = f_itp;
          } else {
            a = b = x_itp;
            fa = fb           = 0;
            this->m_converged = true;
            break;
          }
        }

        this->m_a  = a;
        this->m_b  = b;
        this->m_fa = sign * fa;
        this->m_fb = sign * fb;
        return std::abs(fa) < std::abs(fb) ? a : b;

#undef CMD
      }

    };  // class ITP

  }  // namespace RootFinder

}  // namespace Optimist

#endif  // OPTIMIST_ROOTFINDER_ITP_HH
