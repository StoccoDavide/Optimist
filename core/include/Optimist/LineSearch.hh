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

#ifndef OPTIMIST_LINESEARCH_HH
#define OPTIMIST_LINESEARCH_HH

#include "Optimist.hh"

namespace Optimist {

  /**
   * \brief Namespace for the one-dimensional line searches used by the
   * descent optimizers.
   *
   * \includedoc docs/markdown/LineSearch.md
   *
   * All line searches work on the restriction of the objective to the search
   * direction, \f$ \phi(\alpha) = f(\mathbf{x} + \alpha\mathbf{d}) \f$, and
   * are therefore independent of the vector type. The restriction is provided
   * as a callable with signature
   * \code
   * bool phi(Scalar alpha, Scalar &value, Scalar *derivative);
   * \endcode
   * where \c derivative is \c nullptr when \f$ \phi'(\alpha) \f$ is not
   * needed.
   */
  namespace LineSearch {

    /**
     * \brief Available line search strategies.
     */
    using Method = enum class Method : Integer {
      ARMIJO       = 0, /**< Backtracking with sufficient decrease only. */
      GOLDSTEIN    = 1, /**< Goldstein conditions (no derivatives). */
      STRONG_WOLFE = 2, /**< Strong Wolfe conditions with cubic zoom. */
      MORE_THUENTE = 3  /**< Strong Wolfe conditions, Moré-Thuente search. */
    };

    /**
     * \brief Get the name of a line search strategy.
     */
    inline std::string name(const Method method) {
      switch (method) {
        case Method::ARMIJO:
          return "Armijo";
        case Method::GOLDSTEIN:
          return "Goldstein";
        case Method::STRONG_WOLFE:
          return "StrongWolfe";
        case Method::MORE_THUENTE:
          return "MoreThuente";
        default:
          return "Unknown";
      }
    }

    /**
     * \brief Line search parameters and result.
     * \tparam Scalar Floating-point number type.
     */
    template <typename Scalar>
    struct Settings {
      Scalar c1{1.0e-4};          /**< Sufficient decrease constant. */
      Scalar c2{0.9};             /**< Curvature constant (Wolfe). */
      Scalar goldstein{0.25};     /**< Goldstein constant in (0, 1/2). */
      Scalar shrink{0.5};         /**< Backtracking factor. */
      Scalar alpha_max{1.0e10};   /**< Maximum step. */
      Integer max_iterations{50}; /**< Maximum number of trial steps. */
      Scalar f_slack{
        100 * std::numeric_limits<Scalar>::epsilon()}; /**< Relative slack
            \f$ \epsilon_f \f$ added to the sufficient decrease condition,
            \f$ \phi(\alpha) \leq \phi(0) + c_1\alpha\phi'(0) +
            \epsilon_f|\phi(0)| \f$, so that steps whose decrease is below
            the rounding error of \f$ \phi \f$ are not rejected near the
            minimizer (in the spirit of the approximate Wolfe conditions of
            Hager and Zhang). The factor 100 accounts for cancellation in
            the evaluation of phi. */
      Scalar x_tolerance{
        std::numeric_limits<Scalar>::epsilon()}; /**< Relative tolerance on
                                                    the bracket width
                                                    (Moré-Thuente). */
    };

    /**
     * Right-hand side of the (relaxed) sufficient decrease condition.
     */
    template <typename Scalar>
    Scalar armijo_bound(const Scalar phi_0,
                        const Scalar dphi_0,
                        const Scalar alpha,
                        const Settings<Scalar> &s,
                        const Scalar c) {
      return phi_0 + c * alpha * dphi_0 + s.f_slack * std::abs(phi_0);
    }

    /**
     * Minimizer of the cubic interpolating \f$ (a, \phi_a, \phi'_a) \f$ and
     * \f$ (b, \phi_b, \phi'_b) \f$, safeguarded to the inner part of the
     * interval. Falls back to bisection when the cubic has no minimizer.
     */
    template <typename Scalar>
    Scalar cubic_minimizer(const Scalar a,
                           const Scalar phi_a,
                           const Scalar dphi_a,
                           const Scalar b,
                           const Scalar phi_b,
                           const Scalar dphi_b) {
      const Scalar lo{std::min(a, b)}, hi{std::max(a, b)};
      const Scalar margin{static_cast<Scalar>(0.1) * (hi - lo)};
      const Scalar d1{dphi_a + dphi_b - 3 * (phi_a - phi_b) / (a - b)};
      const Scalar disc{d1 * d1 - dphi_a * dphi_b};
      Scalar x{static_cast<Scalar>(0.5) * (a + b)};
      if (disc >= 0 && std::isfinite(disc)) {
        const Scalar d2{(b > a ? 1 : -1) * std::sqrt(disc)};
        const Scalar den{dphi_b - dphi_a + 2 * d2};
        if (den != 0) {
          const Scalar t{b - (b - a) * (dphi_b + d2 - d1) / den};
          if (std::isfinite(t)) {
            x = t;
          }
        }
      }
      return std::clamp(x, lo + margin, hi - margin);
    }

    /**
     * \brief Backtracking line search enforcing the Armijo condition
     * \f$ \phi(\alpha) \leq \phi(0) + c_1\alpha\phi'(0) \f$, with
     * safeguarded quadratic interpolation of the trial steps.
     * \param[in] phi Restricted objective.
     * \param[in] phi_0 Value \f$ \phi(0) \f$.
     * \param[in] dphi_0 Slope \f$ \phi'(0) < 0 \f$.
     * \param[in,out] alpha Initial trial step on input, accepted step on
     * output.
     * \param[out] phi_alpha Value at the accepted step.
     * \param[in] s Line search settings.
     * \return True if an acceptable step was found.
     */
    template <typename Scalar, typename Phi>
    bool armijo(Phi &&phi,
                const Scalar phi_0,
                const Scalar dphi_0,
                Scalar &alpha,
                Scalar &phi_alpha,
                const Settings<Scalar> &s = Settings<Scalar>()) {
      if (!(dphi_0 < 0)) {
        return false;
      }
      for (Integer i{0}; i < s.max_iterations; ++i) {
        if (!phi(alpha, phi_alpha, nullptr) || !std::isfinite(phi_alpha)) {
          alpha *= s.shrink;
          continue;
        }
        if (phi_alpha <= armijo_bound(phi_0, dphi_0, alpha, s, s.c1)) {
          return true;
        }
        // Minimizer of the quadratic through phi(0), phi'(0) and phi(alpha),
        // safeguarded in [0.1, 0.5] * alpha
        const Scalar den{2 * (phi_alpha - phi_0 - dphi_0 * alpha)};
        Scalar next{den > 0 ? -dphi_0 * alpha * alpha / den : s.shrink * alpha};
        alpha = std::clamp(next,
                           static_cast<Scalar>(0.1) * alpha,
                           static_cast<Scalar>(0.5) * alpha);
      }
      return false;
    }

    /**
     * \brief Line search enforcing the Goldstein conditions
     * \f$ \phi(0) + (1-c)\alpha\phi'(0) \leq \phi(\alpha) \leq \phi(0) +
     * c\alpha\phi'(0) \f$, by expansion and bisection. Only function values
     * are used.
     * \param[in] phi Restricted objective.
     * \param[in] phi_0 Value \f$ \phi(0) \f$.
     * \param[in] dphi_0 Slope \f$ \phi'(0) < 0 \f$.
     * \param[in,out] alpha Initial trial step on input, accepted step on
     * output.
     * \param[out] phi_alpha Value at the accepted step.
     * \param[in] s Line search settings.
     * \return True if an acceptable step was found.
     */
    template <typename Scalar, typename Phi>
    bool goldstein(Phi &&phi,
                   const Scalar phi_0,
                   const Scalar dphi_0,
                   Scalar &alpha,
                   Scalar &phi_alpha,
                   const Settings<Scalar> &s = Settings<Scalar>()) {
      if (!(dphi_0 < 0)) {
        return false;
      }
      const Scalar c{s.goldstein};
      Scalar lo{0}, hi{std::numeric_limits<Scalar>::infinity()};
      for (Integer i{0}; i < s.max_iterations; ++i) {
        const bool ok{phi(alpha, phi_alpha, nullptr) &&
                      std::isfinite(phi_alpha)};
        if (!ok || phi_alpha > armijo_bound(phi_0, dphi_0, alpha, s, c)) {
          hi = alpha;  // Too long
        } else if (phi_alpha < phi_0 + (1 - c) * alpha * dphi_0 -
                                   s.f_slack * std::abs(phi_0)) {
          // Too short (decreases below the rounding level of phi are not
          // significant and are accepted as they are)
          lo = alpha;
        } else {
          return true;
        }
        alpha = std::isfinite(hi) ? static_cast<Scalar>(0.5) * (lo + hi)
                                  : std::min(2 * alpha, s.alpha_max);
      }
      // Accept the last point satisfying the sufficient decrease condition
      if (lo > 0) {
        alpha = lo;
        return phi(alpha, phi_alpha, nullptr) && std::isfinite(phi_alpha);
      }
      return false;
    }

    /**
     * \brief Line search enforcing the strong Wolfe conditions
     * \f$ \phi(\alpha) \leq \phi(0) + c_1\alpha\phi'(0) \f$ and
     * \f$ |\phi'(\alpha)| \leq c_2|\phi'(0)| \f$ (Nocedal and Wright,
     * Algorithms 3.5 and 3.6), with cubic interpolation in the zoom phase.
     * \param[in] phi Restricted objective (with derivative).
     * \param[in] phi_0 Value \f$ \phi(0) \f$.
     * \param[in] dphi_0 Slope \f$ \phi'(0) < 0 \f$.
     * \param[in,out] alpha Initial trial step on input, accepted step on
     * output.
     * \param[out] phi_alpha Value at the accepted step.
     * \param[out] dphi_alpha Slope at the accepted step.
     * \param[in] s Line search settings.
     * \return True if an acceptable step was found.
     */
    template <typename Scalar, typename Phi>
    bool strong_wolfe(Phi &&phi,
                      const Scalar phi_0,
                      const Scalar dphi_0,
                      Scalar &alpha,
                      Scalar &phi_alpha,
                      Scalar &dphi_alpha,
                      const Settings<Scalar> &s = Settings<Scalar>()) {
      if (!(dphi_0 < 0)) {
        return false;
      }
      Integer evaluations{0};
      const Scalar noise{s.f_slack * std::abs(phi_0)};
      auto eval = [&](const Scalar a, Scalar &v, Scalar &d) -> bool {
        ++evaluations;
        return phi(a, v, &d) && std::isfinite(v) && std::isfinite(d);
      };

      // Zoom phase: lo satisfies sufficient decrease and has the lowest value
      // so far, the interval [lo, hi] contains points satisfying the strong
      // Wolfe conditions
      auto zoom = [&](Scalar lo,
                      Scalar phi_lo,
                      Scalar dphi_lo,
                      Scalar hi,
                      Scalar phi_hi,
                      Scalar dphi_hi) -> bool {
        while (evaluations < s.max_iterations) {
          Scalar a{cubic_minimizer(lo, phi_lo, dphi_lo, hi, phi_hi, dphi_hi)};
          Scalar v, d;
          if (!eval(a, v, d)) {
            hi     = a;
            phi_hi = std::numeric_limits<Scalar>::infinity();
            continue;
          }
          // Value comparisons are only meaningful above the rounding level of
          // phi: below it, the derivative decides
          if (v > armijo_bound(phi_0, dphi_0, a, s, s.c1) || v > phi_lo + noise) {
            hi      = a;
            phi_hi  = v;
            dphi_hi = d;
          } else {
            if (std::abs(d) <= -s.c2 * dphi_0) {
              alpha      = a;
              phi_alpha  = v;
              dphi_alpha = d;
              return true;
            }
            if (d * (hi - lo) >= 0) {
              hi      = lo;
              phi_hi  = phi_lo;
              dphi_hi = dphi_lo;
            }
            lo      = a;
            phi_lo  = v;
            dphi_lo = d;
          }
          if (std::abs(hi - lo) <=
              std::numeric_limits<Scalar>::epsilon() * std::max(Scalar(1), lo)) {
            break;
          }
        }
        // Return the best point found if it gives sufficient decrease
        if (lo > 0) {
          alpha = lo;
          return eval(alpha, phi_alpha, dphi_alpha);
        }
        return false;
      };

      // Bracketing phase
      Scalar a_prev{0}, phi_prev{phi_0}, dphi_prev{dphi_0};
      Scalar a{std::min(alpha, s.alpha_max)};
      while (evaluations < s.max_iterations) {
        Scalar v, d;
        if (!eval(a, v, d)) {
          // Undefined function: shrink towards the last good point
          a = a_prev + static_cast<Scalar>(0.5) * (a - a_prev);
          continue;
        }
        if (v > armijo_bound(phi_0, dphi_0, a, s, s.c1) ||
            (a_prev > 0 && v > phi_prev + noise)) {
          return zoom(a_prev, phi_prev, dphi_prev, a, v, d);
        }
        if (std::abs(d) <= -s.c2 * dphi_0) {
          alpha      = a;
          phi_alpha  = v;
          dphi_alpha = d;
          return true;
        }
        if (d >= 0) {
          return zoom(a, v, d, a_prev, phi_prev, dphi_prev);
        }
        a_prev    = a;
        phi_prev  = v;
        dphi_prev = d;
        if (a >= s.alpha_max) {
          break;
        }
        a = std::min(2 * a, s.alpha_max);
      }
      if (a_prev > 0) {
        alpha = a_prev;
        return eval(alpha, phi_alpha, dphi_alpha);
      }
      return false;
    }

    /**
     * \brief Safeguarded step of the Moré-Thuente line search (MINPACK-2
     * \c dcstep). Computes a new trial step from the interval endpoints \f$
     * (s_x, f_x, d_x) \f$, \f$ (s_y, f_y, d_y) \f$ and the current trial \f$
     * (s_p, f_p, d_p) \f$, and updates the interval that contains a step
     * satisfying the strong Wolfe conditions.
     */
    template <typename Scalar>
    void more_thuente_step(Scalar &stx,
                           Scalar &fx,
                           Scalar &dx,
                           Scalar &sty,
                           Scalar &fy,
                           Scalar &dy,
                           Scalar &stp,
                           const Scalar fp,
                           const Scalar dp,
                           bool &brackt,
                           const Scalar stpmin,
                           const Scalar stpmax,
                           const Scalar noise = 0) {
      const Scalar sgnd{dp * (dx / std::abs(dx))};
      Scalar stpf;
      // Function values within the rounding level are considered equal
      const bool higher{fp > fx + noise};
      if (higher) {
        // Case 1: higher function value, the minimum is bracketed
        const Scalar theta{3 * (fx - fp) / (stp - stx) + dx + dp};
        const Scalar s{std::max({std::abs(theta), std::abs(dx), std::abs(dp)})};
        Scalar gamma{s * std::sqrt(std::max(
                             Scalar(0),
                             (theta / s) * (theta / s) - (dx / s) * (dp / s)))};
        if (stp < stx) {
          gamma = -gamma;
        }
        const Scalar p{(gamma - dx) + theta};
        const Scalar q{((gamma - dx) + gamma) + dp};
        const Scalar stpc{stx + (p / q) * (stp - stx)};
        const Scalar stpq{stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) *
                                    (stp - stx)};
        stpf   = std::abs(stpc - stx) < std::abs(stpq - stx)
                     ? stpc
                     : stpc + (stpq - stpc) / 2;
        brackt = true;
      } else if (sgnd < 0) {
        // Case 2: lower function value and derivatives of opposite sign
        const Scalar theta{3 * (fx - fp) / (stp - stx) + dx + dp};
        const Scalar s{std::max({std::abs(theta), std::abs(dx), std::abs(dp)})};
        Scalar gamma{s * std::sqrt(std::max(
                             Scalar(0),
                             (theta / s) * (theta / s) - (dx / s) * (dp / s)))};
        if (stp > stx) {
          gamma = -gamma;
        }
        const Scalar p{(gamma - dp) + theta};
        const Scalar q{((gamma - dp) + gamma) + dx};
        const Scalar stpc{stp + (p / q) * (stx - stp)};
        const Scalar stpq{stp + (dp / (dp - dx)) * (stx - stp)};
        stpf   = std::abs(stpc - stp) > std::abs(stpq - stp) ? stpc : stpq;
        brackt = true;
      } else if (std::abs(dp) < std::abs(dx)) {
        // Case 3: lower function value, same sign, decreasing magnitude
        const Scalar theta{3 * (fx - fp) / (stp - stx) + dx + dp};
        const Scalar s{std::max({std::abs(theta), std::abs(dx), std::abs(dp)})};
        Scalar gamma{s * std::sqrt(std::max(
                             Scalar(0),
                             (theta / s) * (theta / s) - (dx / s) * (dp / s)))};
        if (stp > stx) {
          gamma = -gamma;
        }
        const Scalar p{(gamma - dp) + theta};
        const Scalar q{(gamma + (dx - dp)) + gamma};
        const Scalar r{p / q};
        Scalar stpc;
        if (r < 0 && gamma != 0) {
          stpc = stp + r * (stx - stp);
        } else {
          stpc = stp > stx ? stpmax : stpmin;
        }
        const Scalar stpq{stp + (dp / (dp - dx)) * (stx - stp)};
        if (brackt) {
          stpf = std::abs(stpc - stp) < std::abs(stpq - stp) ? stpc : stpq;
          const Scalar bound{stp + Scalar(0.66) * (sty - stp)};
          stpf = stp > stx ? std::min(bound, stpf) : std::max(bound, stpf);
        } else {
          stpf = std::abs(stpc - stp) > std::abs(stpq - stp) ? stpc : stpq;
          stpf = std::clamp(stpf, stpmin, stpmax);
        }
      } else {
        // Case 4: lower function value, same sign, non-decreasing magnitude
        if (brackt) {
          const Scalar theta{3 * (fp - fy) / (sty - stp) + dy + dp};
          const Scalar s{
            std::max({std::abs(theta), std::abs(dy), std::abs(dp)})};
          Scalar gamma{s * std::sqrt(std::max(Scalar(0),
                                              (theta / s) * (theta / s) -
                                                  (dy / s) * (dp / s)))};
          if (stp > sty) {
            gamma = -gamma;
          }
          const Scalar p{(gamma - dp) + theta};
          const Scalar q{((gamma - dp) + gamma) + dy};
          stpf = stp + (p / q) * (sty - stp);
        } else {
          stpf = stp > stx ? stpmax : stpmin;
        }
      }

      // Update the interval that contains a minimizer
      if (higher) {
        sty = stp;
        fy  = fp;
        dy  = dp;
      } else {
        if (sgnd < 0) {
          sty = stx;
          fy  = fx;
          dy  = dx;
        }
        stx = stp;
        fx  = fp;
        dx  = dp;
      }
      stp = stpf;
    }

    /**
     * \brief Line search of Moré and Thuente (1994), as implemented in
     * MINPACK-2 (\c dcsrch). It finds a step satisfying the strong Wolfe
     * conditions by safeguarded cubic and quadratic interpolation, working on
     * the auxiliary function \f$ \psi(\alpha) = \phi(\alpha) - \phi(0) -
     * c_1\alpha\phi'(0) \f$ until an interval with a suitable step is found.
     * \param[in] phi Restricted objective (with derivative).
     * \param[in] phi_0 Value \f$ \phi(0) \f$.
     * \param[in] dphi_0 Slope \f$ \phi'(0) < 0 \f$.
     * \param[in,out] alpha Initial trial step on input, accepted step on
     * output.
     * \param[out] phi_alpha Value at the accepted step.
     * \param[out] dphi_alpha Slope at the accepted step.
     * \param[in] s Line search settings.
     * \return True if an acceptable step was found.
     */
    template <typename Scalar, typename Phi>
    bool more_thuente(Phi &&phi,
                      const Scalar phi_0,
                      const Scalar dphi_0,
                      Scalar &alpha,
                      Scalar &phi_alpha,
                      Scalar &dphi_alpha,
                      const Settings<Scalar> &s = Settings<Scalar>()) {
      if (!(dphi_0 < 0) || !(alpha > 0)) {
        return false;
      }
      const Scalar stpmin{0}, stpmax{s.alpha_max};
      const Scalar noise{s.f_slack * std::abs(phi_0)};
      const Scalar gtest{s.c1 * dphi_0};
      bool brackt{false}, stage1{true};
      Scalar width{stpmax - stpmin}, width1{2 * width};
      Scalar stx{0}, fx{phi_0}, gx{dphi_0};
      Scalar sty{0}, fy{phi_0}, gy{dphi_0};
      Scalar stp{std::min(alpha, stpmax)};
      Scalar stmin{0}, stmax{stp + 4 * stp};
      Scalar f{phi_0}, g{dphi_0};

      // Best point satisfying the sufficient decrease condition
      Scalar best_stp{0}, best_f{phi_0}, best_g{dphi_0};

      for (Integer i{0}; i < s.max_iterations; ++i) {
        if (!phi(stp, f, &g) || !std::isfinite(f) || !std::isfinite(g)) {
          // Undefined function: move back towards the best endpoint
          stp = stx + (stp - stx) / 2;
          continue;
        }
        const Scalar ftest{armijo_bound(phi_0, dphi_0, stp, s, s.c1)};
        if (f <= ftest && (best_stp == 0 || f < best_f)) {
          best_stp = stp;
          best_f   = f;
          best_g   = g;
        }

        // Convergence: strong Wolfe conditions
        if (f <= ftest && std::abs(g) <= -s.c2 * dphi_0) {
          alpha      = stp;
          phi_alpha  = f;
          dphi_alpha = g;
          return true;
        }

        // Termination by rounding errors or interval collapse
        if ((brackt && (stp <= stmin || stp >= stmax)) ||
            (brackt && stmax - stmin <= s.x_tolerance * stmax) ||
            (stp == stpmax && f <= ftest && g <= gtest) ||
            (stp == stpmin && (f > ftest || g >= gtest))) {
          break;
        }

        // Use the auxiliary function in the first stage, as long as the
        // sufficient decrease condition has not been met with g >= 0
        if (stage1 && f <= ftest && g >= std::min(s.c1, s.c2) * dphi_0) {
          stage1 = false;
        }
        if (stage1 && f <= fx && f > ftest) {
          Scalar fm{f - stp * gtest}, gm{g - gtest};
          Scalar fxm{fx - stx * gtest}, gxm{gx - gtest};
          Scalar fym{fy - sty * gtest}, gym{gy - gtest};
          more_thuente_step(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt,
                            stmin, stmax, noise);
          fx = fxm + stx * gtest;
          fy = fym + sty * gtest;
          gx = gxm + gtest;
          gy = gym + gtest;
        } else {
          more_thuente_step(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin,
                            stmax, noise);
        }

        // Force a sufficient decrease of the interval width
        if (brackt) {
          if (std::abs(sty - stx) >= Scalar(0.66) * width1) {
            stp = stx + (sty - stx) / 2;
          }
          width1 = width;
          width  = std::abs(sty - stx);
          stmin  = std::min(stx, sty);
          stmax  = std::max(stx, sty);
        } else {
          stmin = stp + Scalar(1.1) * (stp - stx);
          stmax = stp + 4 * (stp - stx);
        }
        stp = std::clamp(stp, stpmin, stpmax);
        if ((brackt && (stp <= stmin || stp >= stmax)) ||
            (brackt && stmax - stmin <= s.x_tolerance * stmax)) {
          stp = stx;
        }
      }

      // Fall back to the best point giving sufficient decrease
      if (best_stp > 0) {
        alpha      = best_stp;
        phi_alpha  = best_f;
        dphi_alpha = best_g;
        return true;
      }
      return false;
    }

    /**
     * \brief Run the selected line search.
     *
     * The derivative at the accepted step is returned only by the searches
     * that compute it (strong Wolfe and Moré-Thuente); \c dphi_alpha is left
     * untouched by Armijo and Goldstein.
     * \param[in] method Line search strategy.
     * \param[in] phi Restricted objective.
     * \param[in] phi_0 Value \f$ \phi(0) \f$.
     * \param[in] dphi_0 Slope \f$ \phi'(0) < 0 \f$.
     * \param[in,out] alpha Initial trial step on input, accepted step on
     * output.
     * \param[out] phi_alpha Value at the accepted step.
     * \param[out] dphi_alpha Slope at the accepted step (see above).
     * \param[in] s Line search settings.
     * \return True if an acceptable step was found.
     */
    template <typename Scalar, typename Phi>
    bool search(const Method method,
                Phi &&phi,
                const Scalar phi_0,
                const Scalar dphi_0,
                Scalar &alpha,
                Scalar &phi_alpha,
                Scalar &dphi_alpha,
                const Settings<Scalar> &s) {
      switch (method) {
        case Method::ARMIJO:
          return armijo(phi, phi_0, dphi_0, alpha, phi_alpha, s);
        case Method::GOLDSTEIN:
          return goldstein(phi, phi_0, dphi_0, alpha, phi_alpha, s);
        case Method::MORE_THUENTE:
          return more_thuente(phi, phi_0, dphi_0, alpha, phi_alpha, dphi_alpha,
                              s);
        case Method::STRONG_WOLFE:
        default:
          return strong_wolfe(phi, phi_0, dphi_0, alpha, phi_alpha, dphi_alpha,
                              s);
      }
    }

  }  // namespace LineSearch

}  // namespace Optimist

#endif  // OPTIMIST_LINESEARCH_HH
