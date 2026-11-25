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

#ifndef OPTIMIST_OPTIMIZER_CONJUGATE_GRADIENT_HH
#define OPTIMIST_OPTIMIZER_CONJUGATE_GRADIENT_HH

#include "Optimist/LineSearch.hh"
#include "Optimist/Optimizer.hh"

namespace Optimist {

  namespace Optimizer {

    /**
     * \brief Class container for the nonlinear conjugate gradient method.
     *
     * \includedoc docs/markdown/Optimizer/ConjugateGradient.md
     *
     * \tparam Vector Eigen vector type.
     */
    template <typename Vector>
      requires TypeTrait<Vector>::IsEigen &&
               (!TypeTrait<Vector>::IsFixed || TypeTrait<Vector>::Dimension > 0)
    class ConjugateGradient
        : public Optimizer<Vector, ConjugateGradient<Vector>> {
     public:
      static constexpr bool RequiresFunction{true};
      static constexpr bool RequiresFirstDerivative{true};
      static constexpr bool RequiresSecondDerivative{false};

      using VectorTrait = TypeTrait<Vector>;
      using Scalar      = typename TypeTrait<Vector>::Scalar;
      using typename Optimizer<Vector,
                               ConjugateGradient<Vector>>::FirstDerivative;

      OPTIMIST_BASIC_CONSTANTS(Scalar)

      /**
       * Update formula for the conjugate direction parameter \f$ \beta_k \f$.
       */
      using BetaMethod = enum class BetaMethod : Integer {
        FLETCHER_REEVES    = 0, /**< Fletcher-Reeves. */
        POLAK_RIBIERE      = 1, /**< Polak-Ribière-Polyak. */
        POLAK_RIBIERE_PLUS = 2, /**< Polak-Ribière-Polyak, non-negative. */
        HESTENES_STIEFEL   = 3, /**< Hestenes-Stiefel. */
        CONJUGATE_DESCENT  = 4, /**< Fletcher's conjugate descent. */
        LIU_STOREY         = 5, /**< Liu-Storey. */
        DAI_YUAN           = 6, /**< Dai-Yuan. */
        HAGER_ZHANG        = 7, /**< Hager-Zhang. */
        HAGER_ZHANG_PLUS   = 8  /**< Hager-Zhang, truncated. */
      };

      /**
       * Initial trial step of the line search.
       */
      using AlphaMethod = enum class AlphaMethod : Integer {
        ACCEPTED_STEP    = 0, /**< Last accepted step. */
        BARZILAI_BORWEIN = 1  /**< Barzilai-Borwein step. */
      };

     private:
      BetaMethod m_beta_method{BetaMethod::HAGER_ZHANG}; /**< Beta formula. */
      AlphaMethod m_alpha_method{
        AlphaMethod::ACCEPTED_STEP}; /**< Initial trial step. */
      LineSearch::Method m_line_search{
        LineSearch::Method::STRONG_WOLFE};     /**< Line search strategy. */
      LineSearch::Settings<Scalar> m_settings{
        .c1             = static_cast<Scalar>(1.0e-4),
        .c2             = static_cast<Scalar>(0.1),
        .alpha_max      = 1 / SQRT_EPSILON,
        .max_iterations = 25};                 /**< Line search settings. */
      Scalar m_initial_step{1.0};              /**< Initial trial step. */
      Scalar m_min_step{SQRT_EPSILON};         /**< Minimum trial step. */
      Scalar m_descent_parameter{1.0e-2};      /**< Sufficient descent. */
      Scalar m_truncation_parameter{1.0e-2};   /**< Hager-Zhang+ bound. */
      bool m_powell_restart{true};             /**< Powell restart flag. */
      Scalar m_restart_tolerance{0.9};         /**< Powell restart threshold. */

     public:
      /**
       * Class constructor for the conjugate gradient solver.
       */
      ConjugateGradient() {}

      /**
       * Get the conjugate gradient solver name.
       * \return The conjugate gradient solver name.
       */
      constexpr std::string name_impl() const {
        return "ConjugateGradient";
      }

      /**
       * Get the update formula for \f$ \beta_k \f$.
       */
      BetaMethod beta_method() const {
        return this->m_beta_method;
      }

      /**
       * Set the update formula for \f$ \beta_k \f$.
       */
      void beta_method(const BetaMethod t_beta_method) {
        this->m_beta_method = t_beta_method;
      }

      /** Enable the Fletcher-Reeves formula. */
      void enable_fletcher_reeves_method() {
        this->m_beta_method = BetaMethod::FLETCHER_REEVES;
      }

      /** Enable the Polak-Ribière-Polyak formula. */
      void enable_polak_ribiere_method() {
        this->m_beta_method = BetaMethod::POLAK_RIBIERE;
      }

      /** Enable the non-negative Polak-Ribière-Polyak formula. */
      void enable_polak_ribiere_plus_method() {
        this->m_beta_method = BetaMethod::POLAK_RIBIERE_PLUS;
      }

      /** Enable the Hestenes-Stiefel formula. */
      void enable_hestenes_stiefel_method() {
        this->m_beta_method = BetaMethod::HESTENES_STIEFEL;
      }

      /** Enable Fletcher's conjugate descent formula. */
      void enable_conjugate_descent_method() {
        this->m_beta_method = BetaMethod::CONJUGATE_DESCENT;
      }

      /** Enable the Liu-Storey formula. */
      void enable_liu_storey_method() {
        this->m_beta_method = BetaMethod::LIU_STOREY;
      }

      /** Enable the Dai-Yuan formula. */
      void enable_dai_yuan_method() {
        this->m_beta_method = BetaMethod::DAI_YUAN;
      }

      /** Enable the Hager-Zhang formula. */
      void enable_hager_zhang_method() {
        this->m_beta_method = BetaMethod::HAGER_ZHANG;
      }

      /** Enable the truncated Hager-Zhang formula. */
      void enable_hager_zhang_plus_method() {
        this->m_beta_method = BetaMethod::HAGER_ZHANG_PLUS;
      }

      /**
       * Get the initial trial step strategy.
       */
      AlphaMethod alpha_method() const {
        return this->m_alpha_method;
      }

      /**
       * Set the initial trial step strategy.
       */
      void alpha_method(const AlphaMethod t_alpha_method) {
        this->m_alpha_method = t_alpha_method;
      }

      /** Start each line search from the last accepted step. */
      void enable_accepted_step_alpha_method() {
        this->m_alpha_method = AlphaMethod::ACCEPTED_STEP;
      }

      /** Start each line search from the Barzilai-Borwein step. */
      void enable_barzilai_borwein_alpha_method() {
        this->m_alpha_method = AlphaMethod::BARZILAI_BORWEIN;
      }

      /**
       * Get the line search strategy.
       */
      LineSearch::Method line_search() const {
        return this->m_line_search;
      }

      /**
       * Set the line search strategy. The strong Wolfe conditions (strong
       * Wolfe or Moré-Thuente searches) are required by the convergence
       * theory of most \f$ \beta_k \f$ formulas.
       */
      void line_search(const LineSearch::Method t_line_search) {
        this->m_line_search = t_line_search;
      }

      /**
       * Access the line search settings (by default \f$ c_1 = 10^{-4} \f$ and
       * \f$ c_2 = 0.1 \f$, i.e., a fairly accurate line search).
       */
      LineSearch::Settings<Scalar> &line_search_settings() {
        return this->m_settings;
      }

      /**
       * Get the first trial step.
       */
      Scalar initial_step() const {
        return this->m_initial_step;
      }

      /**
       * Set the first trial step.
       */
      void initial_step(const Scalar t_initial_step) {
        OPTIMIST_ASSERT(std::isfinite(t_initial_step) && t_initial_step > 0,
                        "Optimist::Optimizer::ConjugateGradient::initial_step("
                        "...): invalid input detected.");
        this->m_initial_step = t_initial_step;
      }

      /**
       * Get the minimum trial step.
       */
      Scalar min_step() const {
        return this->m_min_step;
      }

      /**
       * Set the minimum trial step.
       */
      void min_step(const Scalar t_min_step) {
        OPTIMIST_ASSERT(std::isfinite(t_min_step) && t_min_step > 0,
                        "Optimist::Optimizer::ConjugateGradient::min_step(...)"
                        ": invalid input detected.");
        this->m_min_step = t_min_step;
      }

      /**
       * Get the sufficient descent constant \f$ c \f$: a direction is kept
       * only if \f$ \mathbf{g}^\top\mathbf{d} \leq -c\|\mathbf{g}\|^2 \f$.
       */
      Scalar descent_parameter() const {
        return this->m_descent_parameter;
      }

      /**
       * Set the sufficient descent constant \f$ 0 < c < 1 \f$.
       */
      void descent_parameter(const Scalar t_descent_parameter) {
        OPTIMIST_ASSERT(t_descent_parameter > 0 && t_descent_parameter < 1,
                        "Optimist::Optimizer::ConjugateGradient::descent_"
                        "parameter(...): invalid input detected.");
        this->m_descent_parameter = t_descent_parameter;
      }

      /**
       * Get the truncation parameter \f$ \eta \f$ of the Hager-Zhang+ formula.
       */
      Scalar truncation_parameter() const {
        return this->m_truncation_parameter;
      }

      /**
       * Set the truncation parameter \f$ \eta > 0 \f$ of the Hager-Zhang+
       * formula.
       */
      void truncation_parameter(const Scalar t_truncation_parameter) {
        OPTIMIST_ASSERT(t_truncation_parameter > 0,
                        "Optimist::Optimizer::ConjugateGradient::truncation_"
                        "parameter(...): invalid input detected.");
        this->m_truncation_parameter = t_truncation_parameter;
      }

      /**
       * Get the Powell restart flag.
       */
      bool powell_restart() const {
        return this->m_powell_restart;
      }

      /**
       * Set the Powell restart flag: restart from the steepest descent
       * direction when consecutive gradients are far from orthogonal.
       */
      void powell_restart(const bool t_powell_restart) {
        this->m_powell_restart = t_powell_restart;
      }

      /** Enable the Powell restart. */
      void enable_powell_restart() {
        this->m_powell_restart = true;
      }

      /** Disable the Powell restart. */
      void disable_powell_restart() {
        this->m_powell_restart = false;
      }

      /**
       * Get the Powell restart threshold \f$ \nu \f$.
       */
      Scalar restart_tolerance() const {
        return this->m_restart_tolerance;
      }

      /**
       * Set the Powell restart threshold \f$ 0 < \nu < 1 \f$: restart when
       * \f$ \mathbf{g}_{k+1}^\top\mathbf{g}_k \leq -\nu
       * \|\mathbf{g}_{k+1}\|\|\mathbf{g}_k\| \f$.
       */
      void restart_tolerance(const Scalar t_restart_tolerance) {
        OPTIMIST_ASSERT(t_restart_tolerance > 0 && t_restart_tolerance < 1,
                        "Optimist::Optimizer::ConjugateGradient::restart_"
                        "tolerance(...): invalid input detected.");
        this->m_restart_tolerance = t_restart_tolerance;
      }

      /**
       * Solve the minimization problem.
       * \param[in] function Cost function.
       * \param[in] gradient Gradient of the cost function.
       * \param[in] x_ini Initialization point.
       * \param[out] x_sol Solution point.
       * \return The convergence boolean flag.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool solve_impl(FunctionLambda &&function,
                      GradientLambda &&gradient,
                      const Vector &x_ini,
                      Vector &x_sol) {
#define CMD "Optimist::Optimizer::ConjugateGradient::solve_impl(...): "

        this->reset_counters();
        if (this->m_verbose) {
          this->header();
        }

        // Current point (x, f, g), search direction d and trial step
        Vector x(x_ini), x_new(x_ini);
        Scalar f{0}, f_new{0};
        FirstDerivative g, g_new, d;
        OPTIMIST_ASSERT(this->evaluate_function(function, x, f),
                        CMD "function evaluation failed at the initial point.");
        OPTIMIST_ASSERT(this->evaluate_gradient(gradient, x, g),
                        CMD "gradient evaluation failed at the initial point.");
        d = -g;
        Scalar trial{this->clamp_step(this->m_initial_step)}, accepted{0};

        std::string notes{"Initialize"};
        for (this->m_iterations = 1;
             this->m_iterations < this->m_max_iterations;
             ++this->m_iterations) {
          if (this->m_verbose) {
            this->info(g.norm(), notes);
          }
          if (g.norm() < this->m_tolerance) {
            this->m_converged = true;
            break;
          }
          notes = "-";
          if (!this->has_sufficient_descent(g, d)) {
            d     = -g;
            notes = "Restart";
          }

          // Line search. On failure, or if the step stagnates below the
          // rounding level of x, retry once along the steepest descent
          // direction before giving up
          const bool ok{this->line_search(function, gradient, x, f, g, d, trial,
                                          x_new, f_new, g_new, accepted)};
          const bool steepest{(d + g).norm() <= EPSILON * g.norm()};
          const bool stagnated{ok && (x_new - x).norm() <=
                                         EPSILON * (1 + x.norm())};
          if (!ok || stagnated) {
            if (steepest) {
              notes = ok ? "Stagnation" : "Line search failed";
              if (ok) {
                x = x_new;
                f = f_new;
                g = g_new;
                this->m_converged = g.norm() < this->m_tolerance;
              }
              break;
            }
            if (ok) {
              x = x_new;
              f = f_new;
              g = g_new;
            }
            d     = -g;
            trial = this->clamp_step(this->m_initial_step);
            notes = "Restart";
            continue;
          }

          // New conjugate direction
          const FirstDerivative y{g_new - g};
          const Scalar beta{this->compute_beta(g, g_new, d, y)};
          FirstDerivative d_new{-g_new + beta * d};
          if (!std::isfinite(beta) || this->should_restart(g, g_new) ||
              !this->has_sufficient_descent(g_new, d_new)) {
            d_new = -g_new;
            notes = "Restart";
          }
          trial = this->next_trial_step(accepted, d, g, g_new);
          x     = x_new;
          f     = f_new;
          g     = g_new;
          d     = d_new;
        }

        if (this->m_verbose) {
          this->bottom();
        }
        x_sol = x;
        return this->m_converged;

#undef CMD
      }

     private:
      /**
       * Clamp a trial step to \f$ [s_{\min}, \alpha_{\max}] \f$, replacing
       * invalid values with the initial step.
       */
      Scalar clamp_step(Scalar t_step) const {
        if (!std::isfinite(t_step) || t_step <= 0) {
          t_step = this->m_initial_step;
        }
        return std::clamp(t_step, this->m_min_step, this->m_settings.alpha_max);
      }

      /**
       * Initial trial step of the next line search.
       */
      Scalar next_trial_step(const Scalar accepted,
                             const FirstDerivative &d,
                             const FirstDerivative &g,
                             const FirstDerivative &g_new) const {
        if (this->m_alpha_method == AlphaMethod::BARZILAI_BORWEIN) {
          // alpha = s's / s'y along the previous direction, s = accepted * d
          const FirstDerivative y{g_new - g};
          const Scalar dy{d.dot(y)};
          if (dy > EPSILON * d.norm() * y.norm()) {
            return this->clamp_step(accepted * d.squaredNorm() / dy);
          }
        }
        return this->clamp_step(accepted);
      }

      /**
       * Check the sufficient descent condition \f$ \mathbf{g}^\top\mathbf{d}
       * \leq -c\|\mathbf{g}\|^2 \f$.
       */
      bool has_sufficient_descent(const FirstDerivative &g,
                                  const FirstDerivative &d) const {
        const Scalar gd{g.dot(d)};
        return std::isfinite(gd) &&
               gd <= -this->m_descent_parameter * g.squaredNorm();
      }

      /**
       * Powell's restart test on consecutive gradients.
       */
      bool should_restart(const FirstDerivative &g,
                          const FirstDerivative &g_new) const {
        if (!this->m_powell_restart) {
          return false;
        }
        const Scalar n{g.norm() * g_new.norm()};
        return n > 0 && g_new.dot(g) <= -this->m_restart_tolerance * n;
      }

      /**
       * Compute \f$ \beta_k \f$ with the selected formula. Degenerate
       * denominators are detected relative to the scale of the vectors
       * involved, so that the formulas keep working when the gradient is small
       * (i.e., close to the minimizer); in that case \f$ \beta_k = 0 \f$
       * (steepest descent).
       * \param[in] g Previous gradient \f$ \mathbf{g}_k \f$.
       * \param[in] g_new Current gradient \f$ \mathbf{g}_{k+1} \f$.
       * \param[in] d Previous direction \f$ \mathbf{d}_k \f$.
       * \param[in] y Gradient change \f$ \mathbf{y}_k = \mathbf{g}_{k+1} -
       * \mathbf{g}_k \f$.
       * \return The parameter \f$ \beta_k \f$.
       */
      Scalar compute_beta(const FirstDerivative &g,
                          const FirstDerivative &g_new,
                          const FirstDerivative &d,
                          const FirstDerivative &y) const {
        const Scalar gg{g.squaredNorm()};
        const Scalar gg_new{g_new.squaredNorm()};
        const Scalar dg{-d.dot(g)};
        const Scalar dy{d.dot(y)};
        const Scalar gy{g_new.dot(y)};
        const bool ok_gg{gg > 0};
        const bool ok_dg{dg > EPSILON * d.norm() * g.norm()};
        const bool ok_dy{std::abs(dy) > EPSILON * d.norm() * y.norm()};

        // Hager-Zhang: g_new'(y - 2 d |y|^2/(d'y)) / (d'y)
        auto hager_zhang = [&]() -> Scalar {
          return (gy - 2 * y.squaredNorm() / dy * g_new.dot(d)) / dy;
        };

        switch (this->m_beta_method) {
          case BetaMethod::FLETCHER_REEVES:
            return ok_gg ? gg_new / gg : 0;
          case BetaMethod::POLAK_RIBIERE:
            return ok_gg ? gy / gg : 0;
          case BetaMethod::POLAK_RIBIERE_PLUS:
            return ok_gg ? std::max(gy / gg, Scalar(0)) : 0;
          case BetaMethod::HESTENES_STIEFEL:
            return ok_dy ? gy / dy : 0;
          case BetaMethod::CONJUGATE_DESCENT:
            return ok_dg ? gg_new / dg : 0;
          case BetaMethod::LIU_STOREY:
            return ok_dg ? gy / dg : 0;
          case BetaMethod::DAI_YUAN:
            return ok_dy ? gg_new / dy : 0;
          case BetaMethod::HAGER_ZHANG:
            return ok_dy ? hager_zhang() : 0;
          case BetaMethod::HAGER_ZHANG_PLUS: {
            if (!ok_dy) {
              return 0;
            }
            // Truncation eta_k = -1/(|d| min(eta, |g_new|))
            const Scalar den{d.norm() *
                             std::min(this->m_truncation_parameter,
                                      std::sqrt(gg_new))};
            const Scalar bound{den > 0 ? -1 / den : -INFTY};
            return std::max(hager_zhang(), bound);
          }
          default:
            return 0;
        }
      }

      /**
       * Line search along \f$ \mathbf{d} \f$ with the selected strategy.
       * \return True if a step was accepted; in that case \c x_new, \c f_new
       * and \c g_new refer to the accepted point.
       */
      template <typename FunctionLambda, typename GradientLambda>
      bool line_search(FunctionLambda &function,
                       GradientLambda &gradient,
                       const Vector &x,
                       const Scalar f,
                       const FirstDerivative &g,
                       const FirstDerivative &d,
                       const Scalar trial,
                       Vector &x_new,
                       Scalar &f_new,
                       FirstDerivative &g_new,
                       Scalar &accepted) {
#define CMD "Optimist::Optimizer::ConjugateGradient::line_search(...): "

        const Scalar slope{g.dot(d)};
        OPTIMIST_ASSERT(slope < 0, CMD "non-descent search direction.");

        // Restriction of the objective, keeping the last evaluated point
        Scalar last{-1};
        bool g_valid{false};
        this->m_relaxations = 0;
        auto phi = [&](const Scalar a, Scalar &v, Scalar *dv) -> bool {
          ++this->m_relaxations;
          last    = a;
          g_valid = false;
          x_new   = x + a * d;
          if (!this->evaluate_function(function, x_new, v)) {
            return false;
          }
          f_new = v;
          if (dv != nullptr) {
            if (!this->evaluate_gradient(gradient, x_new, g_new)) {
              return false;
            }
            g_valid = true;
            *dv     = g_new.dot(d);
          }
          return true;
        };

        accepted = trial;
        Scalar value{f}, dvalue{slope};
        if (!LineSearch::search(this->m_line_search, phi, f, slope, accepted,
                                value, dvalue, this->m_settings)) {
          return false;
        }

        // Make sure the outputs refer to the accepted step
        if (last != accepted) {
          x_new = x + accepted * d;
          OPTIMIST_ASSERT(this->evaluate_function(function, x_new, f_new),
                          CMD "function evaluation failed.");
          g_valid = false;
        }
        if (!g_valid) {
          OPTIMIST_ASSERT(this->evaluate_gradient(gradient, x_new, g_new),
                          CMD "gradient evaluation failed.");
        }
        return true;

#undef CMD
      }

    };  // class ConjugateGradient

  }  // namespace Optimizer

}  // namespace Optimist

#endif  // OPTIMIST_OPTIMIZER_CONJUGATE_GRADIENT_HH
